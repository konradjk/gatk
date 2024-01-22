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

        // this is a little confusing because we treat the haplotype as a "read"
        final int variantStartInHaplotype = ReadUtils.getReadIndexForReferenceCoordinate(haplotype.getStart(),
                haplotype.getCigar(), vc.getStart()).getLeft();

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
        final List<Integer> result = new ArrayList<>();

        // advance to approximately (lower bound) where the region to encode begins
        int fwdCnt = 0;
        while (!asm.isRightEdge() && asm.getGenomePosition() + padding < vc.getStart()) {
            asm.stepForwardOnGenome();
            fwdCnt++;
            int poo = vc.getStart() - read.getStart() - padding;
        }

        // encode before the variant
        while (!asm.isRightEdge() && asm.getGenomePosition() < vc.getStart()) {
            int encoding = getEncoding(read, refBases, refOffsetOfReadStart, asm);
            beforeVariant.add(encoding);
            asm.stepForwardOnGenome();
        }

        // left edge of read within padding window
        result.addAll(Collections.nCopies(Math.min(padding - beforeVariant.size(), 0), BITSET_PAST_END_OF_READ));

        result.addAll(beforeVariant.subList(Math.max(beforeVariant.size() - padding, 0), beforeVariant.size()));

        // encode (padding + 1) bases at and after the variant
        while (!asm.isRightEdge() && result.size() <= (2 * padding)) {
            int encoding = getEncoding(read, refBases, refOffsetOfReadStart, asm);
            result.add(encoding);
            asm.stepForwardOnGenome();
        }

        // right edge of read past padding window
        result.addAll(Collections.nCopies(Math.max(2*padding + 1 - result.size(), 0), BITSET_PAST_END_OF_READ));
        final StringBuilder builder = new StringBuilder();
        result.forEach(n -> builder.append(String.format("%02X", n)));
        return builder.toString();
    }

    private static int getEncoding(GATKRead read, byte[] refBases, int refOffsetOfReadStart, AlignmentStateMachine asm) {
        int encoding = 0;
        final boolean consumesRead = asm.getCigarOperator().consumesReadBases();
        final boolean consumesRef = asm.getCigarOperator().consumesReferenceBases();

        if (consumesRead) {
            final byte base = read.getBase(asm.getReadOffset());
            final byte qual = read.getBaseQuality(asm.getReadOffset());
            if (base == 'A') {
                encoding += BITSET_A;
            } else if (base == 'C') {
                encoding += BITSET_C;
            } else if (base == 'G') {
                encoding += BITSET_G;
            } else if (base == 'T') {
                encoding += BITSET_T;
            }

            if (consumesRef) {  // match -- check for substitution error
                final byte refBase = refBases[asm.getGenomeOffset() + refOffsetOfReadStart];
                if (refBase != base) {
                    encoding += BITSET_SUBSTITUTION;
                }
            } else {    // insertion
                encoding += BITSET_INSERTION;
            }

            if (qual < 10) {
                encoding += BITSET_VERY_LOW_QUAL;
            } else if (qual < 20) {
                encoding += BITSET_LOW_QUAL;
            }
        } else if (asm.getCigarOperator().consumesReferenceBases()) {   // deletion
            encoding += BITSET_DELETION;
        }
        return encoding;
    }


}
