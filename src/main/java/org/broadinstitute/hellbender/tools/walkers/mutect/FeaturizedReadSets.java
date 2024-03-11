package org.broadinstitute.hellbender.tools.walkers.mutect;


import htsjdk.samtools.CigarElement;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.tools.walkers.annotator.BaseQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReadPosition;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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

    public static final int GOOD_BASE_QUAL = 20;

    public static final int BAD_BASE_QUAL = 10;

    private static final SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);

    private FeaturizedReadSets() { }

    private enum Encoding {
        MATCH('M'),
        PAST_END('E'),
        DELETION('D'),
        INSERTION('I'),
        MISMATCH('X'),
        LOW_QUAL('Q');

        public final Character encoding;

        Encoding( final Character encoding ) {
            this.encoding = encoding;
        }
    }

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


    // reads are encoded as a vector (List<Integer>) and a string summarized the read near the variant (which Permutect later converts to a tensor)
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
        final StringBuilder result = new StringBuilder();

        // advance to approximately (lower bound) where the region to encode begins
        while (!asm.isRightEdge() && (asm.getGenomePosition() + padding < vc.getStart() || asm.isLeftEdge())) {
            asm.stepForwardOnGenome();
        }

        final List<Pair<Encoding, Character>> encodings = new ArrayList<>(2 * padding);

        // we are now within a distance padding and ready to grow the encoding string
        // vc.getStart() +/- padding may be an overestimate and the span we need, since when insertions are encoded we need less reference span
        // thus we may need to trim later
        while (!asm.isRightEdge() && (asm.getGenomePosition() <= vc.getStart() + padding) && (asm.getGenomeOffset() + refOffsetOfReadStart) < refBases.length) {
            final PileupElement pe = asm.makePileupElement();
            final byte base = pe.getBase();
            final byte qual = pe.getQual();
            final byte refBase = refBases[asm.getGenomeOffset() + refOffsetOfReadStart];

            // the stepForwardOnGenome method skips over inserted bases, so we account for those here
            if (pe.isAfterInsertion()) {
                final List<CigarElement> intervening = pe.getBetweenPrevPosition(); // there should only be one element here, the insertion
                final int insertionSize = intervening.get(intervening.size() - 1).getLength();
                final int readOffsetAtInsertionStart = asm.getReadOffset() - insertionSize;

                for (int n = 0; n < insertionSize; n++) {
                    final byte insertedBase = read.getBase(readOffsetAtInsertionStart + n);
                    final byte insertedQual = read.getBaseQuality(readOffsetAtInsertionStart + n);
                    encodings.add(ImmutablePair.of(Encoding.INSERTION, encodeBaseAndQual(insertedBase, insertedQual)));
                }
            }

            // at the event start: trim down to exactly padding encodings and pad with PAST_END if necessary
            if (asm.getGenomePosition() == vc.getStart()) {
                if (encodings.size() < padding) {   // eg E4 if the read start is 4 bases into the padding region
                    result.append(Encoding.PAST_END).append(padding - encodings.size());
                }
                processEncodings(result, encodings.subList(Math.max(encodings.size() - padding, 0), encodings.size()));
                encodings.clear();
            }

            if (pe.isDeletion()) {
                encodings.add(ImmutablePair.of(Encoding.DELETION, ' '));
            } else if (base == refBase) {   // high-qual or low-qual match
                final Encoding encoding = qual >= GOOD_BASE_QUAL ? Encoding.MATCH : Encoding.LOW_QUAL;
                final Character character = qual >= GOOD_BASE_QUAL ? ' ' : encodeBaseAndQual(base, qual);
                encodings.add(ImmutablePair.of(encoding, character));
            } else {    // substitution
                encodings.add(ImmutablePair.of(Encoding.MISMATCH, encodeBaseAndQual(base, qual)));
            }

            asm.stepForwardOnGenome();
        } // done getting encodings, now we must process the encodings at and after the variant

        if (asm.getGenomePosition() == vc.getStart()) {
            // note: padding + 1 comes up because the variant start itself is 1, and after the variant start is padding
            processEncodings(result, encodings.subList(0, Math.min(encodings.size(), padding + 1)));
            if (encodings.size() < padding + 1) {   // eg E4 if the read start is 4 bases into the padding region
                result.append(Encoding.PAST_END).append(padding + 1 - encodings.size());
            }
        }

        return result.toString();
    }

    private static void processEncodings(final StringBuilder result, final List<Pair<Encoding, Character>> outputEncodings) {
        // if we have eg MMMDDMMMM this returns {0,3,5,9}, telling us that runs of identical encoding types are [0,3); [3,5); [5,9)
        final int[] runBounds = IntStream.range(0, outputEncodings.size() + 1)
                .filter(n -> n == 0 || n == outputEncodings.size() || outputEncodings.get(n).getLeft() != outputEncodings.get(n-1).getLeft())
                .toArray();

        // DDD -> D3, MMMM -> M4, QaQc -> Qac, XG -> XG, IAIGIT -> IAGT
        for (int n = 0; n < runBounds.length - 1; n++) {
            final Encoding encoding = outputEncodings.get(runBounds[n]).getLeft();
            final int runLength = runBounds[n+1] - runBounds[n];

            result.append(encoding.encoding);
            if (encoding == Encoding.MATCH || encoding == Encoding.DELETION) {
                result.append(runLength);
            } else {
                new IndexRange(runBounds[n], runBounds[n+1]).forEach(m -> result.append(outputEncodings.get(m).getRight()));
            }
        }
    }


    private static char encodeBaseAndQual(final byte base, final byte qual) {
        if (qual >= GOOD_BASE_QUAL) {
            return (char) base;
        } else if (qual < BAD_BASE_QUAL) {
            return 'N';
        } else  {
            return Character.toLowerCase((char) base);
        }
    }
}
