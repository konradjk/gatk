package org.broadinstitute.hellbender.tools.walkers.mutect;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.tools.walkers.annotator.BaseQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReadPosition;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * For each sample and for each allele a list feature vectors of supporting reads
 * In order to reduce the number of delimiter characters, we flatten featurized reads.  For example, suppose allele 1 has featurized reads
 * [1,2] and [3,4] and allele 2 has featurized reads [5,6] and [7,8], the annotation is
 * 1,2,3,4|5,6,7,8
 */
public class FeaturizedReadSets {
    private static final Logger logger = LogManager.getLogger(FeaturizedReadSets.class);

    public static final int DEFAULT_BASE_QUALITY = 25;

    public static final int PADDING = 10;

    private static final SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
    public static final int BITSET_A = 0;
    public static final int BITSET_C = 1;
    public static final int BITSET_G = 2;
    public static final int BITSET_T = 3;
    public static final int BITSET_DELETION = 4;
    public static final int BITSET_SUBSTITUTION = 8;
    public static final int BITSET_INSERTION = 16;
    public static final int BITSET_LOW_QUAL = 32;
    public static final int BITSET_VERY_LOW_QUAL = 64;
    public static final int BITSET_PAST_END_OF_READ = 128;

    private FeaturizedReadSets() { }

    public static List<List<Pair<List<Integer>, String>>> getReadVectors(final VariantContext vc,
                                                           final Collection<String> samples,
                                                           final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                                           final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods,
                                                           final int refDownsample,
                                                           final int altDownsample,
                                                           final M2ArgumentCollection.Mutect3DatasetMode mutect3DatasetMode) {
        return getReadVectors(vc, samples, likelihoods, haplotypeLikelihoods, refDownsample, altDownsample, Collections.emptyMap(), mutect3DatasetMode);
    }

    // returns Lists (in allele order) of lists of read vectors supporting each allele
    public static List<List<Pair<List<Integer>, String>>> getReadVectors(final VariantContext vc,
                                                           final Collection<String> samples,
                                                           final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                                           final AlleleLikelihoods<Fragment, Haplotype> haplotypeLikelihoods,
                                                           final int refDownsample,
                                                           final int altDownsample,
                                                           final Map<Allele, Integer> altDownsampleMap,
                                                           final M2ArgumentCollection.Mutect3DatasetMode mutect3DatasetMode) {
        final Map<Allele, List<GATKRead>> readsByAllele = likelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> new ArrayList<>()));

        samples.stream().flatMap(s -> likelihoods.bestAllelesBreakingTies(s).stream())
                .filter(ba -> ba.isInformative())
                .forEach(ba -> readsByAllele.get(ba.allele).add(ba.evidence));

        // downsample if necessary
        final Allele refAllele = likelihoods.alleles().stream().filter(Allele::isReference).findFirst().get();
        for (final Allele allele : likelihoods.alleles()) {
            final int downsample = allele.isReference() ? refDownsample : altDownsampleMap.getOrDefault(allele, altDownsample);
            if (readsByAllele.get(allele).size() > downsample) {
                Collections.shuffle(readsByAllele.get(allele));
                readsByAllele.put(allele, readsByAllele.get(allele).subList(0, downsample));
            }
        }

        final Map<GATKRead, Haplotype> bestHaplotypes = new HashMap<>();
        samples.stream().flatMap(s -> haplotypeLikelihoods.bestAllelesBreakingTies(s).stream())
                .forEach(ba -> ba.evidence.getReads().forEach(read -> bestHaplotypes.put(read, ba.allele)));

        return vc.getAlleles().stream()
                .map(allele -> readsByAllele.get(allele).stream().map(read -> featurize(read, vc, bestHaplotypes, mutect3DatasetMode)).collect(Collectors.toList()))
                .collect(Collectors.toList());
    }


    private static Pair<List<Integer>, String> featurize(final GATKRead read, final VariantContext vc,
                                                final Map<GATKRead, Haplotype> bestHaplotypes,
                                                final M2ArgumentCollection.Mutect3DatasetMode mutect3DatasetMode) {
        final List<Integer> result = new ArrayList<>();
        result.add(read.getMappingQuality());
        result.add(BaseQuality.getBaseQuality(read, vc).orElse(DEFAULT_BASE_QUALITY));
        result.add(read.isFirstOfPair() ? 1 : 0);
        result.add(read.isReverseStrand() ? 1 : 0);

        // distances from ends of read
        final int readPosition = ReadPosition.getPosition(read, vc).orElse(0);
        result.add(readPosition);
        result.add(read.getLength() - readPosition);


        result.add(Math.abs(read.getFragmentLength()));

        if (mutect3DatasetMode == M2ArgumentCollection.Mutect3DatasetMode.ILLUMINA) {
            // distances from ends of fragment
            final int fragmentStart = Math.min(read.getMateStart(), read.getUnclippedStart());
            final int fragmentEnd = fragmentStart + Math.abs(read.getFragmentLength());
            result.add(vc.getStart() - fragmentStart);
            result.add(fragmentEnd - vc.getEnd());
        }

        // Ultima-specific read tags
        if (mutect3DatasetMode == M2ArgumentCollection.Mutect3DatasetMode.ULTIMA) {
            result.add(read.getAttributeAsInteger("si"));   // si is an integer on the order of 100s or 1000s
            result.add((int) (1000*read.getAttributeAsFloat("rq")));    // rq is a float from 0 and 1, so we multiply by 1000 and round
        }

        // mismatches versus best haplotype
        final Haplotype haplotype = bestHaplotypes.get(read);

        // TODO: fix this
        // I have no idea why this edge case occurs in Ultima data and maybe sometimes in Illumina
        if (!bestHaplotypes.containsKey(read)) {
            logger.warn(String.format("Best haplotypes don't contain read %s at position %s:%d.", read.getName(),
                    vc.getContig(), vc.getStart()));
            result.add(3);
            result.add(2);
            return null;
        }

        final SmithWatermanAlignment readToHaplotypeAlignment = aligner.align(haplotype.getBases(), read.getBases(), SmithWatermanAlignmentConstants.ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, SWOverhangStrategy.SOFTCLIP);
        final int readStartInHaplotype = readToHaplotypeAlignment.getAlignmentOffset();

        final GATKRead copy = read.copy();
        copy.setCigar(readToHaplotypeAlignment.getCigar());
        final int mismatchCount = AlignmentUtils.getMismatchCount(copy, haplotype.getBases(), readStartInHaplotype).numMismatches;
        result.add(mismatchCount);

        final long indelsVsBestHaplotype = readToHaplotypeAlignment.getCigar().getCigarElements().stream().filter(el -> el.getOperator().isIndel()).count();
        result.add((int) indelsVsBestHaplotype);

        Utils.validate(result.size() == mutect3DatasetMode.getNumReadFeatures(), "Wrong number of features");

        final String encodingAroundVariant = encodeVicinityOfVariant(read, vc, haplotype.getBases(),
                readStartInHaplotype, PADDING);

        return Pair.of(result, encodingAroundVariant);
    }


    /**
     * We encode reads near a variant as follows.  First, we imagine lining up the read's bases with those of the best
     * haplotype like so:
     *
     *  read:        ACGTNNNACGGGGACGTXXXX
     *  haplotype:   ACGTTTTACGNNNACGTAACC
     *
     * this means the read has a TTTT -> T deletion, a G -> GGGG insertion, and ends before the final AACC.
     *
     * Then for some number of bases before and after the variant we encode each position.  If we use 5 bases, for example,
     * and the second C above is the variant position then we would encode the part
     *
     *  read:        TNNNACGGGGA
     *  haplotype:   TTTTACGNNNA
     *
     * Note that this means neither 5 reference bases nor 5 read bases, but five CIGAR bases in which we count elements
     * that consume read OR reference bases.
     *
     * The encoding is like a bit set.  First we have a part from 0 to 7 for the read bases: 0 for A, 1 for C, 2 for G,
     * 3 for T, 4 for deletion, 5 for past beginning/end of read.  Then add 16 for an insertion, 8 for a substitution error,
     *
     *
     * @param padding number of bases used strictly before (ref coordinate is less than refCoord) given ref position, and
     *                coincident OR after it.
     * @param read
     * @param refBases  bases of ref sequence
     * @param refOffsetOfReadStart index of read soft start within ref byte[] array
     */
    public static String encodeVicinityOfVariant(final GATKRead read, final VariantContext vc, final byte[] refBases, final int refOffsetOfReadStart, final int padding) {
        final AlignmentStateMachine asm = new AlignmentStateMachine(read);  // starts off the left edge, not on the read

        // this list may be built bigger than it needs to be if there are insertions just before the variant
        final List<Integer> beforeVariant = new ArrayList<>();
        final List<Integer> afterVariant = new ArrayList<>();

        // advance to approximately (lower bound) where the region to encode begins
        while (!asm.isRightEdge() && (asm.getGenomePosition() + padding < vc.getStart() || asm.isLeftEdge())) {
            asm.stepForwardOnGenome();
        }

        while (!asm.isRightEdge() && afterVariant.size() <= padding && (asm.getGenomeOffset() + refOffsetOfReadStart) < refBases.length) {
            final List<Integer> listToGrow = (asm.getGenomePosition() < vc.getStart()) ? beforeVariant : afterVariant;

            final PileupElement pe = asm.makePileupElement();
            // the stepForwardOnGenome method skips over inserted bases, so we check for them before advancing
            if (pe.isBeforeInsertion()) {
                final int insertionLength = pe.getLengthOfImmediatelyFollowingIndel();
                for (int n = 0; n < insertionLength; n++) {

                    final byte base = read.getBase(asm.getReadOffset() + n + 1);
                    final byte qual = read.getBaseQuality(asm.getReadOffset() + n + 1);
                    listToGrow.add(bitsetEncodingOfBase(base) + bitsetEncodingOfBaseQual(qual) + BITSET_INSERTION);
                }
            }
            listToGrow.add(getEncoding(read, refBases, refOffsetOfReadStart, asm));
            asm.stepForwardOnGenome();
        }

        final List<Integer> result = new ArrayList<>();
        // left edge of read within padding window
        result.addAll(Collections.nCopies(Math.max(padding - beforeVariant.size(), 0), BITSET_PAST_END_OF_READ));
        result.addAll(beforeVariant.subList(Math.max(beforeVariant.size() - padding, 0), beforeVariant.size()));
        result.addAll(afterVariant.subList(0, Math.min(afterVariant.size(), padding + 1)));
        result.addAll(Collections.nCopies(Math.max(padding + 1 - afterVariant.size(), 0), BITSET_PAST_END_OF_READ));

        Utils.validate(result.size() == 2 * padding + 1, () -> "wrong result size");
        final StringBuilder builder = new StringBuilder();
        result.forEach(n -> builder.append(String.format("%02X", n)));  // 2 hexadecimal digits with leading zero if relevant
        return builder.toString();
    }

    private static int getEncoding(GATKRead read, byte[] refBases, int refOffsetOfReadStart, AlignmentStateMachine asm) {
        if (asm.getCigarOperator().consumesReadBases()) {   // match/substitution
            final byte base = read.getBase(asm.getReadOffset());
            final byte qual = read.getBaseQuality(asm.getReadOffset());
            final byte refBase = refBases[asm.getGenomeOffset() + refOffsetOfReadStart];

            return bitsetEncodingOfBase(base) + bitsetEncodingOfBaseQual(qual) + (refBase == base ? 0 : BITSET_SUBSTITUTION);
        } else {   // deletion
            return BITSET_DELETION;
        }
    }

    private static int bitsetEncodingOfBaseQual(byte qual) {
        if (qual > 20) {
            return 0;
        } else if (qual > 10) {
            return BITSET_LOW_QUAL;
        } else {
            return BITSET_VERY_LOW_QUAL;
        }
    }

    private static int bitsetEncodingOfBase(byte base) {
        if (base == 'A') {
            return BITSET_A;
        } else if (base == 'C') {
            return BITSET_C;
        } else if (base == 'G') {
            return BITSET_G;
        } else if (base == 'T') {
            return BITSET_T;
        } else {
            return 0;
        }
    }
}
