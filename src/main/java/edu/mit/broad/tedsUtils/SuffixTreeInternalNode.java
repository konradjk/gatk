/*
 * $Id: SuffixTreeInternalNode.java 73459 2008-10-10 16:10:03Z tsharpe $
 * WHITEHEAD INSTITUTE
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2002 by the
 * Whitehead Institute for Biomedical Research.  All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support
 * whatsoever.  The Whitehead Institute can not be responsible for its
 * use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;

/**
 * An non-leaf node in a suffix tree.
 *
 * @author tsharpe
 * @version $Revision: 48897 $
 */
public class SuffixTreeInternalNode
    extends SuffixTreeLeafNode
{
    public SuffixTreeInternalNode( CharSequence str, int labelStart, int labelEnd )
    {
        super(str, labelStart);
        mLabelEnd = labelEnd;
    }

    public SuffixTreeInternalNode( CharSequence str, int labelStart, int labelEnd, Object strings )
    {
        super(str, labelStart, strings);
        mLabelEnd = labelEnd;
    }

    @Override
    public int getLabelEnd()
    {
        return mLabelEnd;
    }

    @Override
    public int getLabelLength()
    {
        return mLabelEnd - mLabelStart;
    }

    @Override
    public SuffixTreeNode getChild( char chr )
    {
        SuffixTreeNode result = null;
        chr = Character.toUpperCase(chr);

        if ( chr == 'A' )
        {
            result = mALink;
        }
        else if ( chr == 'C' )
        {
            result = mCLink;
        }
        else if ( chr == 'G' )
        {
            result = mGLink;
        }
        else if ( chr == 'T' )
        {
            result = mTLink;
        }

        return result;
    }

    @Override
    public void setChild( char chr, SuffixTreeNode child )
    {
        chr = Character.toUpperCase(chr);
        if ( chr == 'A' )
        {
            mALink = child;
        }
        else if ( chr == 'C' )
        {
            mCLink = child;
        }
        else if ( chr == 'G' )
        {
            mGLink = child;
        }
        else if ( chr == 'T' )
        {
            mTLink = child;
        }
        else
        {
            throw new IllegalArgumentException("Illegal character: " + chr);
        }
    }

    @Override
    public SuffixTreeNode getSuffixLink()
    {
        return mSuffixLink;
    }

    @Override
    public void setSuffixLink( SuffixTreeNode suffixLink )
    {
        mSuffixLink = suffixLink;
    }

    @Override
    public boolean isLeaf()
    {
        return false;
    }

    @Override
    public boolean isRoot()
    {
        return mLabelEnd == 0;
    }

    private int mLabelEnd;
    private SuffixTreeNode mALink;
    private SuffixTreeNode mCLink;
    private SuffixTreeNode mGLink;
    private SuffixTreeNode mTLink;
    private SuffixTreeNode mSuffixLink;
}
