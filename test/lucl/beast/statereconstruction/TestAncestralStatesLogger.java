/**
 * 
 */
package lucl.beast.statereconstruction;

import org.junit.Assert;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import beast.evolution.tree.Node;
import beast.util.TreeParser;

/**
 * @author kaipingga
 *
 */
class TestAncestralStatesLogger {
	/**
	 * Test method for
	 * {@link lucl.beast.statereconstruction.AncestralStatesLogger#initAndValidate()}.
	 */
	@Test
	void testInitAndValidate() {
		Assert.fail("Not yet implemented");
	}

	/**
	 * Test method for
	 * {@link lucl.beast.statereconstruction.AncestralStatesLogger#init(java.io.PrintStream)}.
	 */
	@Test
	void testInit() {
		Assert.fail("Not yet implemented");
	}

	/**
	 * Test method for
	 * {@link lucl.beast.statereconstruction.AncestralStatesLogger#log(int, java.io.PrintStream)}.
	 */
	@Test
	void testLog() {
		Assert.fail("Not yet implemented");
	}

	/**
	 * Test method for
	 * {@link lucl.beast.statereconstruction.AncestralStatesLogger#commonAncestors(java.util.List)}.
	 */
	@Test
	void testCommonAncestors() {
		boolean isLabeled = true;
		TreeParser parser = new TreeParser("((A:1,B:1)C:2,((D:1,E:1)F:1,G:2)H:1):1", false, false, isLabeled, 1);
		Node root = parser.getRoot();
		Node c = root.getChild(0);
		Node h = root.getChild(1);
		Node f = h.getChild(0);
		Node d = f.getChild(0);
		Node e = f.getChild(1);
		Node g = h.getChild(1);

		List<Node> focusNodes = new ArrayList<Node>(4);
		List<Node> commonAncestors;

		focusNodes.clear();
		focusNodes.add(c);
		commonAncestors = AncestralStatesLogger.commonAncestors(focusNodes);
		Assert.assertEquals(c, commonAncestors.get(commonAncestors.size() - 1));

		focusNodes.clear();
		focusNodes.add(c);
		focusNodes.add(h);
		commonAncestors = AncestralStatesLogger.commonAncestors(focusNodes);
		Assert.assertEquals(root, commonAncestors.get(commonAncestors.size() - 1));

		focusNodes.clear();
		focusNodes.add(e);
		focusNodes.add(g);
		commonAncestors = AncestralStatesLogger.commonAncestors(focusNodes);
		Assert.assertEquals(h, commonAncestors.get(commonAncestors.size() - 1));
		Assert.assertFalse(commonAncestors.contains(c));

		focusNodes.clear();
		focusNodes.add(d);
		focusNodes.add(e);
		focusNodes.add(g);
		commonAncestors = AncestralStatesLogger.commonAncestors(focusNodes);
		Assert.assertEquals(h, commonAncestors.get(commonAncestors.size() - 1));
		Assert.assertEquals(root, commonAncestors.get(0));
		
		focusNodes.clear();
		focusNodes.add(e);
		focusNodes.add(h);
		commonAncestors = AncestralStatesLogger.commonAncestors(focusNodes);
		Assert.assertEquals(h, commonAncestors.get(commonAncestors.size() - 1));
		Assert.assertFalse(commonAncestors.contains(c));
	}

}
