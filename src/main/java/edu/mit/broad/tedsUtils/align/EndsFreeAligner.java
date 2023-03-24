/**
 * $Id: EndsFreeAligner.java 73459 2008-10-10 16:10:03Z tsharpe $
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
import edu.mit.broad.tedsUtils.align.Alignment.BlockType;

/**
 * Aligns sequence X against sequence Y without penalizing starting or ending gaps in either sequence.
 * In other words, the optimal prefix of either sequence is matched against the optimal suffix of the other one. 
 * For best memory efficiency sequence X ought to be the shorter of the two.
 *
 * @author tsharpe
 * @version $Revision$
 */
public class EndsFreeAligner
    extends Aligner
{
    public EndsFreeAligner( CharSequence seqX, CharSequence seqY )
    {
        super(seqX,seqY,new Scorer(),false);
    }

    public EndsFreeAligner( CharSequence seqX, CharSequence seqY, Scorer scorer )
    {
        super(seqX,seqY,scorer,false);
    }

    public EndsFreeAligner( CharSequence seqX, CharSequence seqY, Scorer scorer, boolean lateGapping )
    {
        super(seqX,seqY,scorer,lateGapping);
    }

    @Override
    public Score getScore()
    {
        return getScore(new Score(),0.F,0.F,0.F,0.F);
    }

    @Override
    public Alignment getAlignment()
    {
        Score score = new Score();
        LinkedList<Block> blocks = getAlignment(score,0.F,0.F,0.F,0.F);
        if ( blocks.getFirst().getType() != BlockType.PAIRWISE )
        {
            blocks.removeFirst();
        }
        return new Alignment(score.getScore(),blocks,new CharSubSequence(mSeqX,0,mSeqX.length()),new CharSubSequence(mSeqY,0,mSeqY.length()));
    }

    @Override
    public Aligner clone( CharSequence seqX, CharSequence seqY )
    {
        return new EndsFreeAligner(seqX,seqY,getScorer(),isLateGapping());
    }
}
