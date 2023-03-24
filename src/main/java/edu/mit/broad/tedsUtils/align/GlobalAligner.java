/**
 * $Id: GlobalAligner.java 73459 2008-10-10 16:10:03Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package edu.mit.broad.tedsUtils.align;

import java.util.LinkedList;

import edu.mit.broad.tedsUtils.CharSubSequence;
import edu.mit.broad.tedsUtils.align.Alignment.Block;

/**
 * Sequence X, in its entirety, is optimally aligned to sequence Y, in its entirety.
 *
 * @author tsharpe
 * @version $Revision$
 */
public class GlobalAligner
    extends Aligner
{
    public GlobalAligner( CharSequence seqX, CharSequence seqY )
    {
        super(seqX,seqY,new Scorer(),false);
    }

    public GlobalAligner( CharSequence seqX, CharSequence seqY, Scorer scorer )
    {
        super(seqX,seqY,scorer,false);
    }

    public GlobalAligner( CharSequence seqX, CharSequence seqY, Scorer scorer, boolean lateGapping )
    {
        super(seqX,seqY,scorer,lateGapping);
    }

    @Override
    public Score getScore()
    {
        AffineGapScorer xGap = mScorer.getXGapScorer();
        AffineGapScorer yGap = mScorer.getYGapScorer();
        return getScore(new GlobalScore(),xGap.getOpenScore(),xGap.getExtendScore(),yGap.getOpenScore(),yGap.getExtendScore());
    }

    @Override
    public Alignment getAlignment()
    {
        AffineGapScorer xGap = mScorer.getXGapScorer();
        AffineGapScorer yGap = mScorer.getYGapScorer();
        Score score = new GlobalScore();
        LinkedList<Block> blocks = getAlignment(score,xGap.getOpenScore(),xGap.getExtendScore(),yGap.getOpenScore(),yGap.getExtendScore());
        return new Alignment(score.getScore(),blocks,new CharSubSequence(mSeqX,0,mSeqX.length()),new CharSubSequence(mSeqY,0,mSeqY.length()));
    }

    @Override
    public Aligner clone( CharSequence seqX, CharSequence seqY )
    {
        return new GlobalAligner(seqX,seqY,getScorer(),isLateGapping());
    }

    static class GlobalScore
        extends Score
    {
        @Override
        protected void checkRowEnd( float score, int xIdx, int yIdx )
        {
            // row ends aren't interesting -- nothing to do
        }

        @Override
        protected void checkLastRow( float[] scores, int yIdx )
        {
            int xIdx = scores.length - 1; // check only the final cell
            checkScore(scores[xIdx],xIdx,yIdx);
        }
    }
}
