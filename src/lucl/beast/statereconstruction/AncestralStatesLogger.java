/**
 * 
 */
package lucl.beast.statereconstruction;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.UserDataType;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * @author kaipingga
 *
 */

@Description("Logs internal states sampled from the distribution at the MRCAs of multiple sets of taxa")
public class AncestralStatesLogger extends TreeLikelihood implements Loggable {
	public Input<List<TaxonSet>> taxonsetInput = new Input<List<TaxonSet>>("taxonset",
			"set of taxa defining a clade. The MRCA node of the clade is logged", new ArrayList<TaxonSet>());
	public Input<String> valueInput = new Input<>("value",
			"space delimited set of labels, one for each site in the alignment. Used as site label in the log file.");
	public Input<Boolean> logParentInput = new Input<>("logParent",
			"flag to indicate the parent value should be logged", false);
	private HashMap<Node, Integer[]> samples = new HashMap<Node, Integer[]>();

	@Override
	public void initAndValidate() {
		// ensure we do not use BEAGLE
		boolean forceJava = Boolean.valueOf(System.getProperty("java.only"));
		System.setProperty("java.only", "true");
		super.initAndValidate();
		System.setProperty("java.only", "" + forceJava);

		String values = valueInput.get();
		if (values != null && values.trim().length() > 0) {
			// use values as labels
			values = values.trim().replaceAll("\\s+", "\t");
			String[] strs = values.split("\t");
			if (strs.length != dataInput.get().getSiteCount()) {
				throw new IllegalArgumentException("Number of labels (" + strs.length
						+ ") does not match amountof data (" + dataInput.get().getSiteCount() + ") " + values);
			}
		}

		List<String> taxaNames = dataInput.get().getTaxaNames();
		for (final TaxonSet set : taxonsetInput.get()) {
			Set<String> taxaSet = new LinkedHashSet<>();
			taxaSet.clear();
			for (final String sTaxon : set.asStringList()) {
				final int iTaxon = taxaNames.indexOf(sTaxon);
				if (iTaxon < 0) {
					throw new IllegalArgumentException("Cannot find taxon " + sTaxon + " in data");
				}
				if (taxaSet.contains(sTaxon)) {
					throw new IllegalArgumentException(
							"Taxon " + sTaxon + " is defined multiple times, while they should be unique");
				}
				if (logParentInput.get() && taxaNames.size() == set.getTaxonCount()) {
					throw new RuntimeException(
							"Cannot log parent of the root; either choose a different clade, or set logParent flag to false");
				}
				taxaSet.add(sTaxon);
			}
		}
	}

	private static String formattedState(Integer[] state, DataType dataType) {
		// The idea of this function is taken from the beast-classic ASTL
		StringBuffer sb = new StringBuffer();
		if (dataType instanceof UserDataType) {
			for (int i : state) {
				sb.append(dataType.getCode(i));
				sb.append("\t");
			}
		} else {
			for (int i : state) {
				sb.append(dataType.getChar(i));
				sb.append("\t");
			}
		}
		return sb.toString();
	}

	@Override
	public void init(PrintStream out) {
		for (TaxonSet taxa : taxonsetInput.get()) {
			String idstring = taxa.getID();
			String values = valueInput.get();
			if (values != null && values.trim().length() > 0) {
				// use values as labels
				values = values.trim().replaceAll("\\s+", "\t" + idstring + ".");
				out.append(values);
				out.append("\t");
			} else {
				int siteCount = dataInput.get().getSiteCount();
				for (int i = 0; i < siteCount; i++) {
					out.append("idstring." + i + "\t");
				}
			}
		}
	}

	@Override
	public void log(int nSample, PrintStream out) {
		// force fresh recalculation of likelihood at this stage
		Arrays.fill(m_branchLengths, 0);
		calculateLogP();

		samples.clear();
		for (TaxonSet taxonset : taxonsetInput.get()) {
			// determine the MRCA node we are going to log
			List<String> taxa = taxonset.asStringList();
			List<Node> leaves = treeInput.get().getRoot().getAllLeafNodes();
			for (Node leaf : leaves) {
				if (!taxa.contains(leaf.getID())) {
					leaves.remove(leaf);
				}
			}
			List<Node> common = commonAncestors(leaves);

			// sample states
			Integer[] sampled;
			if (logParentInput.get()) {
				sampled = sample(common.get(common.size() - 2));
			} else {
				sampled = sample(common.get(common.size() - 1));
			}
			// generate output: convert output according to data type
			out.append(formattedState(sampled, dataInput.get().getDataType()));
		}
	}

	/**
	 * traverse to the root then, sample root values, and propagate back to the MRCA
	 * along the path that goes between root and MRCA
	 * 
	 * @return sample
	 */
	private Integer[] sample(Node node) {
		Integer[] sample = samples.get(node);
		if (sample != null) {
			return sample;
		}
		int siteCount = dataInput.get().getSiteCount();
		int stateCount = dataInput.get().getMaxStateCount();
		sample = new Integer[siteCount];

		if (node.isRoot()) {
			if (beagle != null) {
				throw new RuntimeException("BEAGLE is not supported yet");
				// m_fRootPartials = beagle.m_fRootPartials;
			}

			double[] p = new double[stateCount];

			for (int i = 0; i < sample.length; i++) {
				int offset = stateCount * dataInput.get().getPatternIndex(i);
				for (int j = 0; j < stateCount; j++) {
					p[j] = m_fRootPartials[offset + j];
				}
				sample[i] = Randomizer.randomChoicePDF(p);
			}

		} else {
			Integer[] parentSample = sample(node.getParent());

			double[] p = new double[stateCount];
			double[] partials = new double[dataInput.get().getPatternCount() * stateCount
					* m_siteModel.getCategoryCount()];

			if (m_siteModel.getCategoryCount() != 1) {
				throw new RuntimeException("Gamma rate heterogeneity or proportion invariant is not supported yet");
			}
			if (beagle != null) {
				throw new RuntimeException("BEAGLE is not supported yet");
				// beagle.beagle.getPartials(arg0, arg1, arg2);
				// getTransitionMatrix(nodeNum, probabilities);
			} else {
				likelihoodCore.getNodeMatrix(node.getNr(), 0, probabilities);
			}

			if (node.isLeaf()) {
				if (!m_useAmbiguities.get()) {
					// leaf node values come mainly from the states.
					// only ambiguous sites are sampled

					int[] states = new int[dataInput.get().getPatternCount() * m_siteModel.getCategoryCount()];
					if (beagle != null) {
						throw new RuntimeException("BEAGLE is not supported yet");
						// beagle.beagle.getPartials(arg0, arg1, arg2);
						// getTransitionMatrix(nodeNum, probabilities);
					} else {
						likelihoodCore.getNodeStates(node.getNr(), states);
					}

					for (int j = 0; j < sample.length; j++) {
						int childIndex = dataInput.get().getPatternIndex(j);
						if (states[childIndex] >= 0 && states[childIndex] < stateCount) {
							// copy state, if it is not ambiguous
							sample[j] = states[childIndex];
						} else {
							sample[j] = -1;
						}
					}
				} else {
					// useAmbiguities == true
					// sample conditioned on child partials
					likelihoodCore.getNodePartials(node.getNr(), partials);

					// sample using transition matrix and parent states
					for (int j = 0; j < sample.length; j++) {
						int parentIndex = parentSample[j] * stateCount;
						int childIndex = dataInput.get().getPatternIndex(j) * stateCount;

						for (int i = 0; i < stateCount; i++) {
							p[i] = partials[childIndex + i] * probabilities[parentIndex + i];
						}

						sample[j] = Randomizer.randomChoicePDF(p);
					}
				}
			} else {
				likelihoodCore.getNodePartials(node.getNr(), partials);

				// sample using transition matrix and parent states
				for (int j = 0; j < sample.length; j++) {
					int parentIndex = parentSample[j] * stateCount;
					int childIndex = dataInput.get().getPatternIndex(j) * stateCount;

					for (int i = 0; i < stateCount; i++) {
						p[i] = partials[childIndex + i] * probabilities[parentIndex + i];
					}

					sample[j] = Randomizer.randomChoicePDF(p);
				}
			}
		}
		samples.put(node, sample);
		return sample;
	}

	@Override
	public void close(PrintStream out) {
		// nothing to do
	}

	public static List<Node> commonAncestors(final List<Node> taxa) {
		if (taxa.size() == 0) {
			throw new IllegalArgumentException("No taxa given");
		} else {
			Node here = taxa.get(0);
			List<Node> routeToRoot = new ArrayList<>();
			routeToRoot.add(here);
			while (!here.isRoot()) {
				here = here.getParent();
				routeToRoot.add(0, here);
			}

			if (taxa.size() == 1) {
				return routeToRoot;
			}

			int i = 0;
			List<Node> commonRoute = new ArrayList<>(routeToRoot.size());
			for (Node other : commonAncestors(taxa.subList(1, taxa.size()))) {
				if (routeToRoot.size() < i) {
					return routeToRoot;
				}
				if (other == routeToRoot.get(i)) {
					commonRoute.add(other);
				} else {
					return commonRoute;
				}
				i += 1;
			}
			return commonRoute;
		}
	}

}
