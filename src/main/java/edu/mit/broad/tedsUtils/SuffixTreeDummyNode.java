/**
 * $Id: SuffixTreeDummyNode.java 73459 2008-10-10 16:10:03Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package edu.mit.broad.tedsUtils;

import java.util.Iterator;

/**
 * A fake node.
 * There are certain places in the SuffixTree algorithm where we refer to a node
 * that appears to the compiler as if it may never have been initialized.  Nonetheless,
 * the algorithm will, in fact, always set it.  To check that assertion, as well as to
 * avoid compiler warnings, we initialize the variable to this dummy node which will
 * explode if it's ever used.
 *
 * @author tsharpe
 * @version $Revision$
 */
public class SuffixTreeDummyNode
    implements SuffixTreeNode
{
    public CharSequence getString()
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public char charAt( int idx )
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public int getLabelStart()
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public int getLabelEnd()
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public int getLabelLength()
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public CharSequence getLabel()
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public SuffixTreeNode split( int idx )
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public SuffixTreeNode internalize( int labLen, CharSequence str, int labStart )
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public SuffixTreeNode getChild( char chr )
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public void setChild( char chr, SuffixTreeNode child )
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public SuffixTreeNode getSuffixLink()
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public void setSuffixLink( SuffixTreeNode suffixLink )
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public Iterator<CharSequence> getStringIterator()
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public void addString( CharSequence str )
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public boolean ubiquitousMatch()
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public boolean isRoot()
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    public boolean isLeaf()
    {
        throw new IllegalStateException("This dummy node was never supposed to be referenced.");
    }

    private SuffixTreeDummyNode()
    {
        // nothing to do
    }

    public static final SuffixTreeNode SINGLETON = new SuffixTreeDummyNode();
}
