package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import com.google.api.client.util.Lists;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.apache.commons.lang.math.IntRange;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.LongHomopolymerHaplotypeCollapsingEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.ChainPruner;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;

import java.nio.ByteBuffer;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class AlignmentAugmentedAssembler {

    // TODO: magic constant, make adjustable
    private static final int MAX_BRANCH_VERTICES = 100;

    private final int kmerSize;

    private final byte minBaseQualityToUseInAssembly;

    public AlignmentAugmentedAssembler(final int kmerSize, final byte minBaseQualityToUseInAssembly) {
        this.kmerSize = kmerSize;
        this.minBaseQualityToUseInAssembly = minBaseQualityToUseInAssembly;
    }

    public AssemblyResultSet runLocalAssembly(final AssemblyRegion assemblyRegion, final Haplotype refHaplotype,
                                              final byte[] fullReferenceWithPadding, final SimpleInterval refLoc,
                                              final Iterable<GATKRead> reads,
                                              final ChainPruner<MultiDeBruijnVertex, MultiSampleEdge> pruner,
                                              final SmithWatermanAligner aligner,
                                              final LongHomopolymerHaplotypeCollapsingEngine haplotypeCollapsing,
                                              final SWParameters haplotypeToReferenceSWParameters) {
        // make and prune a simple de Bruijn graph in order to know which kmers are worth tracking
        final PlainDeBruijnGraph initialGraph = new PlainDeBruijnGraph(kmerSize, minBaseQualityToUseInAssembly);
        initialGraph.addSequence("ref", refHaplotype.getBases(), 1, true);
        reads.forEach(read -> initialGraph.addRead(read, null));
        initialGraph.buildGraphIfNecessary();
        pruner.pruneLowWeightChains(initialGraph);

        final Set<Kmer> kmersAfterPruning = initialGraph.vertexSet().stream()
                .map(vertex -> new Kmer(vertex.getSequence()))
                .collect(Collectors.toSet());

        final List<Pair<Kmer, IntRange>> rangedKmers = Utils.stream(reads).flatMap(read -> kmerizeRead(read).stream())
                .filter(pair -> pair.getLeft() != null).toList();

        final List<Pair<Kmer, Integer>> unambiguousKmers = rangedKmers.stream()
                .filter(pair -> kmersAfterPruning.contains(pair.getLeft()))
                .filter(pair -> Math.abs(pair.getRight().getMaximumInteger() - pair.getRight().getMinimumInteger()) < kmerSize)
                .map(pair -> Pair.of(pair.getLeft(), (pair.getRight().getMinimumInteger() + pair.getRight().getMaximumInteger())/2))
                .collect(Collectors.toList());
        unambiguousKmers.addAll(kmerizeReference(refHaplotype));

        final VertexManager vertexManager = new VertexManager(kmerSize, unambiguousKmers);

        final AugmentedKmerGraph graph = makeAugmentedKmerGraph(refHaplotype, reads, vertexManager);

        // edge case: it is possible that the vertex manager may contain a vertex not contained in the graph
        // when this occurs the kmer sequence of the bad vertex exists in the graph at a different position
        final List<AugmentedVertex> bad = vertexManager.allVertices().filter(v -> !graph.containsVertex(v)).toList();
        bad.forEach(v -> vertexManager.removeVertex(v));

        // kind of hacky way to deal with homopolymers / STRs causing identical nearby kmers to collapse onto
        // a single vertex
        if (graph.hasCycles()) {
            return null;
        }

        final List<AugmentedVertex> branchVertices = graph.vertexSet().stream()
                .filter(v -> graph.outDegreeOf(v) > 1).toList();

        // decision vertices = anything after a branch or a source (wince which source to begin with is a decision)
        final List<AugmentedVertex> decisionVertices = Stream.concat(graph.getSources().stream(),
                branchVertices.stream().flatMap(bv -> graph.outgoingVerticesOf(bv).stream()))
                .distinct() // needed because a decision vertex may have multiplo parents
                .toList();

        final Map<AugmentedVertex, Integer> decisionVertexIndexMap = IntStream.range(0, decisionVertices.size()).boxed()
                .collect(Collectors.toMap(decisionVertices::get, n->n));

        final List<Integer> sourceVertexIndices = graph.getSources().stream()
                .map(decisionVertexIndexMap::get).toList();
        
        final Map<GATKRead, double[]> readFeatureMap = Utils.stream(reads)
                .collect(Collectors.toMap(read -> read, read -> featurizeRead(vertexManager, graph, decisionVertices,
                        decisionVertexIndexMap, sourceVertexIndices, read)));

        // the above map with the +1 elements zeroed out
        final Map<GATKRead, double[]> readNegativeFeatureMap = Utils.stream(reads)
                .collect(Collectors.toMap(read -> read, read -> MathUtils.applyToArray(readFeatureMap.get(read), x -> Math.min(x,0))));

        List<double[]> haplotypeCentroids = new ArrayList<>();
        haplotypeCentroids.add(new double[decisionVertices.size()]);

        final Set<GATKRead> readsToSeedNewCluster = new HashSet<>();

        final Random rand = new Random();
        for (int iteration = 0; iteration < 10; iteration++) {
            final List<ArrayList<GATKRead>> readsByHaplotype = haplotypeCentroids.stream().map(hc -> new ArrayList<GATKRead>()).collect(Collectors.toList());
            // assign reads to haplotype via cosine similarity, breaking ties randomly
            for (final GATKRead read : reads) {
                if (!readsToSeedNewCluster.contains(read)) {
                    final double[] features = readFeatureMap.get(read);
                    final double[] haplotypeCosineSimilarities =
                            haplotypeCentroids.stream().mapToDouble(centroid -> MathUtils.dotProduct(centroid, features)).toArray();
                    int randomBestHaplotypeIndex = randomMaxIndex(rand, haplotypeCosineSimilarities);
                    readsByHaplotype.get(randomBestHaplotypeIndex).add(read);
                }
            }

            if (!readsToSeedNewCluster.isEmpty()) {
                readsByHaplotype.add(new ArrayList<>(readsToSeedNewCluster));
                readsToSeedNewCluster.clear();
            }

            final List<double[]> haplotypeMeans = new ArrayList<>();
            for (final List<GATKRead> readsInCluster : readsByHaplotype) {
                final int count = readsInCluster.size();
                final double[] mean = new double[decisionVertices.size()];
                if (count > 0) {
                    readsInCluster.forEach(read -> MathUtils.addToArrayInPlace(mean, readFeatureMap.get(read)));
                    MathUtils.applyToArrayInPlace(mean, x -> x / count);
                }
                haplotypeMeans.add(mean);
            }

            // now traverse the graph to convert double[] haplotype means into int[] features encoding decisions
            haplotypeCentroids.clear();
            for (final double[] haplotypeMean : haplotypeMeans) {
                final double[] centroid = new double[decisionVertices.size()];

                // choose starting source vertex
                final int startVertexIndex = randomMaxIndex(rand, sourceVertexIndices.stream().mapToDouble(n -> haplotypeMean[n]).toArray());
                sourceVertexIndices.forEach(n -> centroid[n] = -1);
                centroid[startVertexIndex] = 1;

                AugmentedVertex vertex = decisionVertices.get(startVertexIndex);
                while (!graph.isSink(vertex)) {

                    if (graph.outDegreeOf(vertex) == 1) {    // zip ahead until the next branch or sink
                        vertex = graph.getEdgeTarget(graph.outgoingEdgeOf(vertex));
                    } else if (graph.outDegreeOf(vertex) > 1) {    // it's a branch
                        final List<AugmentedVertex> children = new ArrayList<>(graph.outgoingVerticesOf(vertex));

                        final double[] childrensMeans = children.stream()
                                .mapToDouble(c -> haplotypeMean[decisionVertexIndexMap.get(c)]).toArray();

                        final AugmentedVertex bestChild = children.get(randomMaxIndex(rand, childrensMeans));

                        children.forEach(c -> centroid[decisionVertexIndexMap.get(c)] = -1);
                        centroid[decisionVertexIndexMap.get(bestChild)] = 1;

                        vertex = bestChild;
                    }
                }   // done traversing graph and featurizing the haplotype centroid
                haplotypeCentroids.add(centroid);
            }

            // decide which haplotype is the best candidate for splitting
            double maxOfMaxDisagreements = 0;
            int clusterToSplit = 0;
            int indexToSplitBy = 0;
            for (int n = 0; n < haplotypeCentroids.size(); n++) {
                final double[] centroid = haplotypeCentroids.get(n);
                final double[] centroidPositive = MathUtils.applyToArray(centroid, x -> Math.max(x,0));
                final List<GATKRead> readsInCluster = readsByHaplotype.get(n);

                final double[] totalNegative = new double[centroid.length];
                readsInCluster.forEach(read -> MathUtils.addToArrayInPlace(totalNegative, readNegativeFeatureMap.get(read)));

                // the entries of this array are the counts of reads in the cluster that reject (-1 entry) a decision
                // vertex chosen (+1 entry) by the centroid haplotype
                final double[] disagreement = MathUtils.applyToArrayInPlace(MathArrays.ebeMultiply(centroidPositive, totalNegative), x -> -x);
                final double maxDisagreement = MathUtils.arrayMax(disagreement);

                if (maxDisagreement >= maxOfMaxDisagreements) {
                    maxOfMaxDisagreements = maxDisagreement;
                    clusterToSplit = n;
                    indexToSplitBy = MathUtils.maxElementIndex(disagreement);
                }
            }

            // seed new cluster with reads belonging to the cluster to be split that disagree at the most "controversial"
            // decision vertex
            int finalIndexToSplitBy = indexToSplitBy;
            readsByHaplotype.get(clusterToSplit).stream()
                    .filter(read -> readFeatureMap.get(read)[finalIndexToSplitBy] == -1)
                    .forEach(readsToSeedNewCluster::add);

            final Set<Integer> emptyAndRedundantClusters = IntStream.range(0, readsByHaplotype.size())
                    .filter(n -> readsByHaplotype.get(n).isEmpty()).boxed().collect(Collectors.toSet());

            final Set<double[]> uniqueCentroids = new HashSet<>();
            for (int n = 0; n < haplotypeCentroids.size(); n++) {
                final double[] centroid = haplotypeCentroids.get(n);
                if (!uniqueCentroids.add(centroid)) {
                    emptyAndRedundantClusters.add(n);
                }
            }

            // remove empty and redundant clusters
            haplotypeCentroids = IntStream.range(0, haplotypeCentroids.size())
                    .filter(n -> !emptyAndRedundantClusters.contains(n))
                    .mapToObj(haplotypeCentroids::get)
                    .collect(Collectors.toList());

            // done with iteration -- we have assigned reads to centroids and recomputed centroids
        }

        // at this point the centroids are arrays of 1, -1, and 0, and at each branch in the haplotype (including when choosing a source)
        // exactly one of the outgoing vertices has an entry of +1, with the rest having -1
        final List<Haplotype> haplotypes = new ArrayList<>(haplotypeCentroids.size());
        for (final double[] centroid : haplotypeCentroids) {
            final List<AugmentedVertex> vertices = new ArrayList<>();

            final int startVertexIndex = sourceVertexIndices.stream()
                    .max(Comparator.comparingDouble(n -> centroid[n])).get();
            AugmentedVertex vertex = decisionVertices.get(startVertexIndex);
            vertices.add(vertex);

            // append a byte (base) after advancing each vertex
            while (!graph.isSink(vertex)) {
                if (graph.outDegreeOf(vertex) == 1) {
                    vertex = graph.getEdgeTarget(graph.outgoingEdgeOf(vertex));
                } else if (graph.outDegreeOf(vertex) > 1) {    // it's a branch; choose the best child
                    vertex = graph.outgoingVerticesOf(vertex).stream()
                            .max(Comparator.comparingDouble(c -> centroid[decisionVertexIndexMap.get(c)])).get();
                }
                vertices.add(vertex);
            }
            // the source vertex adds a full kmer; all others add a base
            final ByteBuffer baseBuffer = ByteBuffer.allocate(vertices.size() + kmerSize - 1);
            baseBuffer.put(vertices.get(0).getSequence());
            IntStream.range(1, vertices.size()).forEach(n -> baseBuffer.put(vertices.get(n).getSuffix()));

            // TODO: need to check whether it matches the reference! reference
            // TODO: probably do that by comparing the centroid byte[]

            // TODO: Lines 358- of ReadThreadingAssembler deal with some edge cases that may come up. . .
            final SmithWatermanAlignment alignment = aligner.align(refHaplotype.getBases(), baseBuffer.array(),
                    haplotypeToReferenceSWParameters, SWOverhangStrategy.SOFTCLIP);

            // TODO: Fix this edge case: if a haplotype has a dangling end that extends PAST the reference haplotype
            // TODO: ie it starts before and/or ends after the ref haplotype this alignment will yield leading or
            // TODO: trailing softclips, which will cause the haplotype trimming later to throw an error.
            // TODO: therefore we need to trim those extra bases and the softclip here

            final Cigar cigar = alignment.getCigar();
            final int leadingSoftclips = cigar.isLeftClipped() ? cigar.getFirstCigarElement().getLength() : 0;
            final int trailingSoftclips = cigar.isRightClipped() ? cigar.getLastCigarElement().getLength() : 0;
            final boolean clippingNeeded = (leadingSoftclips + trailingSoftclips > 0);

            final byte[] trimmedArray =  !clippingNeeded ? baseBuffer.array() :
                    Arrays.copyOfRange(baseBuffer.array(), leadingSoftclips, baseBuffer.array().length - trailingSoftclips);

            final SmithWatermanAlignment trimmedAlignment = !clippingNeeded ? alignment :
                    aligner.align(refHaplotype.getBases(), trimmedArray, haplotypeToReferenceSWParameters, SWOverhangStrategy.SOFTCLIP);

            // after aligning haplotype to ref (remember, we have done nothing for dangling heads/tails, so the haplotype
            // may start after / end before the reference) we add any leading and trailing ref bases and pad the CIGAR
            // with corresponding M sequences.
            final int leadingRefPadding = trimmedAlignment.getAlignmentOffset();
            final int trailingRefPadding = refHaplotype.getBases().length - trimmedAlignment.getCigar().getReferenceLength() - leadingRefPadding;

            final Cigar paddedCigar = new CigarBuilder()
                    .add(new CigarElement(leadingRefPadding, CigarOperator.M))
                    .addAll(trimmedAlignment.getCigar().getCigarElements())
                    .add(new CigarElement(trailingRefPadding, CigarOperator.M))
                    .make();

            final ByteBuffer paddedBaseBuffer = ByteBuffer.allocate(trimmedArray.length + leadingRefPadding + trailingRefPadding);
            paddedBaseBuffer.put(refHaplotype.getBases(), 0, leadingRefPadding);
            paddedBaseBuffer.put(trimmedArray);
            paddedBaseBuffer.put(refHaplotype.getBases(), refHaplotype.getBases().length - trailingRefPadding, trailingRefPadding);

            final byte[] paddedBases = paddedBaseBuffer.array();
            final boolean matchesRef = paddedCigar.numCigarElements() == 1 && Arrays.equals(paddedBases, refHaplotype.getBases());

            final Haplotype haplotype = new Haplotype(paddedBases, matchesRef);
            haplotype.setKmerSize(kmerSize);
            haplotype.setCigar(paddedCigar);
            haplotype.setAlignmentStartHapwrtRef(refHaplotype.getAlignmentStartHapwrtRef());
            haplotype.setGenomeLocation(refHaplotype.getGenomeLocation());
            haplotypes.add(haplotype);
        }


        final AssemblyResultSet resultSet = new AssemblyResultSet();
        resultSet.setRegionForGenotyping(assemblyRegion);
        resultSet.setFullReferenceWithPadding(fullReferenceWithPadding);
        resultSet.setPaddedReferenceLoc(refLoc);
        final SimpleInterval activeRegionExtendedLocation = assemblyRegion.getPaddedSpan();
        refHaplotype.setGenomeLocation(activeRegionExtendedLocation);
        resultSet.add(refHaplotype);
        resultSet.setHaplotypeCollapsingEngine(haplotypeCollapsing);
        haplotypes.forEach(resultSet::add);

        return resultSet;



    }

    private int randomMaxIndex(final Random rand, final double[] array) {
        final double maxValue = MathUtils.arrayMax(array);
        final int[] maxIndices = IntStream.range(0, array.length).filter(n -> array[n] == maxValue).toArray();
        return maxIndices[rand.nextInt(maxIndices.length)];
    }

    private double[] featurizeRead(final VertexManager vertexManager, final AugmentedKmerGraph graph, final List<AugmentedVertex> decisionVertices,
                                final Map<AugmentedVertex, Integer> decisionVertexIndexMap, final List<Integer> sourceVertexIndices, GATKRead read) {
        AugmentedVertex lastVertexInRead = null;
        final double[] features = new double[decisionVertices.size()];

        for (final Pair<Kmer, IntRange> pair : kmerizeRead(read)) {
            final Kmer kmer = pair.getLeft();
            final IntRange positionRange = pair.getRight();
            final Optional<AugmentedVertex> maybeVertex = kmer == null ? Optional.empty() : vertexManager.getVertex(kmer, positionRange);

            if (maybeVertex.isPresent()) {
                final AugmentedVertex vertex = maybeVertex.get();
                // if last vertex was null due to being at the read start or bad bases disqualifying a kmer in
                // the middle of the read, traverse backwards until the last decision vertex
                if (lastVertexInRead == null) {
                    AugmentedVertex backwardVertex = vertex;

                    // traverse backwards as long as it is possible (in-degree > 0) and unambiguous (in-degree < 2)
                    // until we reach 1) a source vertex 2) a non-source decision vertex or 3) a vertex with in-degree > 1
                    // Cases 1 and 2 give us an unambiguous entry in the feature vector; case 3 does not.
                    while (graph.inDegreeOf(backwardVertex) == 1 && !decisionVertexIndexMap.containsKey(backwardVertex)) {
                        backwardVertex = graph.getEdgeSource(graph.incomingEdgeOf(backwardVertex));
                    }

                    if (graph.isSource(backwardVertex)) {
                        sourceVertexIndices.forEach(n -> features[n] = -1);
                        features[decisionVertexIndexMap.get(backwardVertex)] = 1;
                    } else if (graph.inDegreeOf(backwardVertex) == 1) { // in rare edge case of two parents, it's ambiguous and we leave the feature at zero
                        final AugmentedVertex parent = graph.getEdgeSource(graph.incomingEdgeOf(backwardVertex));
                        graph.outgoingVerticesOf(parent).stream()
                                .filter(decisionVertexIndexMap::containsKey)
                                .map(decisionVertexIndexMap::get)
                                .forEach(n -> features[n] = -1);
                        features[decisionVertexIndexMap.get(backwardVertex)] = 1;
                    }
                } else if (decisionVertexIndexMap.containsKey(vertex)) {  // decision vertex
                    // set sibling features to -1 and this feature to +1
                    graph.outgoingVerticesOf(lastVertexInRead).forEach(v -> features[decisionVertexIndexMap.get(v)] = -1);
                    features[decisionVertexIndexMap.get(vertex)] = 1;
                }
            }

            lastVertexInRead = maybeVertex.orElse(null);
        }
        return features;
    }

    private AugmentedKmerGraph makeAugmentedKmerGraph(final Haplotype refHaplotype, Iterable<GATKRead> reads, VertexManager vertexManager) {
        final AugmentedKmerGraph graph = new AugmentedKmerGraph(kmerSize);
        vertexManager.allVertices().forEach(graph::addVertex);

        final List<Pair<Kmer, Integer>> refKmers = kmerizeReference(refHaplotype);
        for (int n = 0; n < refKmers.size() - 1; n++) {
            final Kmer kmer1 = refKmers.get(n).getLeft();
            final Kmer kmer2 = refKmers.get(n+1).getLeft();

            final AugmentedVertex vertex1 = vertexManager.getVertex(kmer1, refKmers.get(n).getRight()).get();
            final AugmentedVertex vertex2 = vertexManager.getVertex(kmer2, refKmers.get(n+1).getRight()).get();
            graph.addEdge(vertex1, vertex2);
        }


        // TODO: code duplication between ref and reads!!!
        for (final GATKRead read : reads) {
            final List<Pair<Kmer, IntRange>> kmers = kmerizeRead(read);
            for (int n = 0; n < kmers.size() - 1; n++) {
                final Kmer kmer1 = kmers.get(n).getLeft();
                final Kmer kmer2 = kmers.get(n+1).getLeft();

                if (kmer1 != null && kmer2 != null) {   // null if unusable kmers due to 'N' or low-qual bases
                    final Optional<AugmentedVertex> vertex1 = vertexManager.getVertex(kmer1, kmers.get(n).getRight());
                    final Optional<AugmentedVertex> vertex2 = vertexManager.getVertex(kmer2, kmers.get(n+1).getRight());

                    if (vertex1.isPresent() && vertex2.isPresent()) {
                        graph.addEdge(vertex1.get(), vertex2.get());
                    }
                }
            }
        }

        graph.removeSingletonOrphanVertices();
        return graph;
    }

    // extracts kmers and their alignment positions -- soft clips are assigned a range with one ambiguous end eg
    // [-infinity, 10] or [10, infinity]
    private List<Pair<Kmer, IntRange>> kmerizeRead(final GATKRead read) {
        final byte[] sequence = read.getBases();
        final byte[] qualities = read.getBaseQualities();
        int numKmers = sequence.length - kmerSize + 1;

        if (numKmers < 1) {
            return List.of();
        }
        final List<Pair<Kmer, IntRange>> result = new ArrayList<>(numKmers);

        final int start = read.hasAttribute(ReadUtils.ORIGINAL_SOFTCLIP_START_TAG) ?
                read.getAttributeAsInteger(ReadUtils.ORIGINAL_SOFTCLIP_START_TAG) : read.getStart();
        final int end = read.hasAttribute(ReadUtils.ORIGINAL_SOFTCLIP_END_TAG) ?
                read.getAttributeAsInteger(ReadUtils.ORIGINAL_SOFTCLIP_END_TAG) : read.getEnd();

        int lastUnusableReadOffset = -1;   // last kmer-disqualifying base within the current kmer span
        for (int n = Math.min(kmerSize, sequence.length) - 1; n >= 0; n--) {
            if (!baseIsUsableForAssembly(sequence[n], qualities[n])) {
                lastUnusableReadOffset = n;
                break;
            }
        }

        int readOffset = 0;
        int refPosition = read.getSoftStart();
        for (final CigarElement cigarElement : read.getCigarElements()) {
            final CigarOperator op = cigarElement.getOperator();
            final int length = cigarElement.getLength();

            if (!op.consumesReadBases()) {  // deletion: skip ahead on reference but don't make kmers
                refPosition += op.consumesReferenceBases() ? length : 0;
            } else {    // make kmers
                for (int n = 0; n < length; n++) {
                    if (readOffset + kmerSize > sequence.length) {  // kmer would go off end of read
                        break;
                    } else {
                        final int endOfKmer = readOffset + kmerSize - 1;
                        if (!baseIsUsableForAssembly(sequence[endOfKmer], qualities[endOfKmer])) {
                            lastUnusableReadOffset = endOfKmer;
                        }

                        if (lastUnusableReadOffset < readOffset) {  // make a kmer
                            final Kmer kmer = new Kmer(sequence, readOffset, kmerSize);
                            final int rangeMin = refPosition < start ? Integer.MIN_VALUE : refPosition;
                            // TODO: verify that this is > and not >=
                            final int rangeMax = refPosition > end ? Integer.MAX_VALUE : refPosition;
                            result.add(Pair.of(kmer, new IntRange(rangeMin, rangeMax)));
                        } else {
                            result.add(Pair.of(null, new IntRange(refPosition)));   // unusable kmer
                        }

                        refPosition += op.consumesReferenceBases() ? 1 : 0;
                        readOffset++;
                    }
                }
            }
        }
        return result;
    }

    private List<Pair<Kmer, Integer>> kmerizeReference(final Haplotype refHaplotype) {
        final byte[] sequence = refHaplotype.getBases();
        final int refStart = refHaplotype.getStart();
        return IntStream.range(0, sequence.length - kmerSize + 1)
                .mapToObj(offset -> Pair.of(new Kmer(sequence, offset, kmerSize), refStart + offset ))
                .toList();
    }

    private static class VertexManager {
        private final Map<Kmer, PositionClusterer> kmerMap = new HashMap<>();

        public VertexManager(final int tolerance, final Iterable<Pair<Kmer, Integer>> kmersAndPositions) {
            for (final Pair<Kmer, Integer> pair : kmersAndPositions) {
                kmerMap.computeIfAbsent(pair.getLeft(), kmer -> new PositionClusterer(tolerance, kmer)).add(pair.getRight());
            }
            kmerMap.values().forEach(clusterer -> clusterer.finish());
        }

        public void removeVertex(final AugmentedVertex v) {
            final Kmer kmer = new Kmer(v.getSequence());
            final PositionClusterer clusterer = kmerMap.get(kmer);
            if (clusterer != null) {
                clusterer.removeVertex(v);
            }
        }

        public Optional<AugmentedVertex> getVertex(final Kmer kmer, final IntRange range) {
            final int rangeMin = range.getMinimumInteger();
            final int rangeMax = range.getMaximumInteger();

            // TODO: magic conventions!
            if (rangeMin < 0) {
                return getVertexFromMaxPosition(kmer, rangeMax);
            } else if (rangeMax == Integer.MAX_VALUE) {
                return getVertexFromMinPosition(kmer, rangeMin);
            } else {
                return getVertex(kmer, (rangeMin + rangeMax)/2);
            }
        }

        public Optional<AugmentedVertex> getVertex(final Kmer kmer, final int position) {
            final PositionClusterer clusterer = kmerMap.get(kmer);
            return clusterer == null ? Optional.empty() : clusterer.getVertex(position);
        }

        public Optional<AugmentedVertex> getVertexFromMaxPosition(final Kmer kmer, final int maxPosition) {
            final PositionClusterer clusterer = kmerMap.get(kmer);
            return clusterer == null ? Optional.empty() : clusterer.getVertexFromMaxPosition(maxPosition);
        }

        public Optional<AugmentedVertex> getVertexFromMinPosition(final Kmer kmer, final int minPosition) {
            final PositionClusterer clusterer = kmerMap.get(kmer);
            return clusterer == null ? Optional.empty() : clusterer.getVertexFromMinPosition(minPosition);
        }

        public Stream<AugmentedVertex> allVertices() {
            return kmerMap.values().stream().flatMap(clusterer -> clusterer.getVertices().stream());
        }

    }

    private static class PositionClusterer {
        private final Kmer kmer;

        // this will always get sorted from greatest to least count whenever thing get too wrong
        private final List<Pair<Integer, MutableInt>> clustersAndCounts = new ArrayList<>(1);

        private final int tolerance;

        private boolean finalized;

        private final List<AugmentedVertex> vertices = new ArrayList<>(1);


        public PositionClusterer(final int tolerance, final Kmer kmer) {
            this.tolerance = tolerance;
            this.kmer = kmer;
        }

        public List<AugmentedVertex> getVertices() { return vertices; }

        public void removeVertex(final AugmentedVertex v) {
            if (Arrays.equals(v.getSequence(), kmer.bases())) {
                final OptionalInt matchingIndex = IntStream.range(0, vertices.size())
                        .filter(n -> vertices.get(n).equals(v))
                        .findFirst();
                matchingIndex.ifPresent(n -> vertices.remove(n));
            }
        }

        public void finish() {
            if (finalized) {
                return;
            }
            finalized = true;
            sort();
            clustersAndCounts.forEach(pair -> vertices.add(new AugmentedVertex(kmer.bases(), pair.getLeft())));
        }

        public Optional<AugmentedVertex> getVertex(final int position) {
            if (vertices.size() == 1) {
                return Optional.of(vertices.get(0));
            } else {
                for (int n = 0; n < vertices.size(); n++) {
                    if (inCluster(position, vertices.get(n).getPosition())) {
                        return Optional.of(vertices.get(n));
                    }
                }
            }
            return Optional.empty();
        }

        public Optional<AugmentedVertex> getVertexFromMaxPosition(final int maxPosition) {
            int matchIndex = -1;
            for (int n = 0; n < vertices.size(); n++) {
                if (vertices.get(n).getPosition() >= maxPosition) {
                    if (matchIndex != -1) {
                        return Optional.empty();    // multiple matches, ambiguous soft clip
                    }
                    matchIndex = n;
                }
            }

            return matchIndex != -1 ? Optional.of(vertices.get(matchIndex)) : Optional.empty();
        }

        public Optional<AugmentedVertex> getVertexFromMinPosition(final int minPosition) {
            int matchIndex = -1;
            for (int n = 0; n < vertices.size(); n++) {
                if (vertices.get(n).getPosition() >= minPosition) {
                    if (matchIndex != -1) {
                        return Optional.empty();    // multiple matches, ambiguous soft clip
                    }
                    matchIndex = n;
                }
            }

            return matchIndex != -1 ? Optional.of(vertices.get(matchIndex)) : Optional.empty();
        }



        public void add(final int position) {
            Utils.validate(!finalized, "already finalized");
            boolean newCluster = true;
            for (int n = 0; n < clustersAndCounts.size(); n++) {
                final MutableInt count = clustersAndCounts.get(n).getRight();
                if (inCluster(position, clustersAndCounts.get(n).getLeft())) {
                    count.increment();
                    if (n > 0 && count.getValue() > 2 * clustersAndCounts.get(0).getRight().getValue()) {
                        sort();
                    }
                    newCluster = false;
                    break;
                }
            }

            if (newCluster) {
                clustersAndCounts.add(Pair.of(position, new MutableInt(1)));
            }
        }

        private void sort() {
            clustersAndCounts.sort(Comparator.comparingInt(pair -> -pair.getRight().getValue()));
        }

        private boolean inCluster(final int position, final int cluster) {
            return Math.abs(position - cluster) < tolerance;
        }
    }

    protected boolean baseIsUsableForAssembly(final byte base, final byte qual) {
        return base != BaseUtils.Base.N.base && qual >= minBaseQualityToUseInAssembly;
    }
}
