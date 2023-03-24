/**
 * $Id: AffineGapScorer.java 57447 2008-01-28 17:41:24Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package edu.mit.broad.tedsUtils.align;

/**
 * Affine model for gap scoring.
 *
 * @author tsharpe
 * @version $Revision$
 */
public class AffineGapScorer
{
    /**
     * Make one.  Both args should be non-positive numbers.
     */
    public AffineGapScorer( float openScore, float extendScore )
    {
        if ( openScore > 0.F || extendScore > 0.F )
        {
            throw new IllegalArgumentException("Affine gap scores must be penalties (<= 0.)");
        }
        mOpenScore = openScore;
        mExtendScore = extendScore;
    }

    /**
     * Penalty for opening a gap.
     */
    public float getOpenScore()
    {
        return mOpenScore;
    }

    /**
     * Penalty for extending a gap.
     */
    public float getExtendScore()
    {
        return mExtendScore;
    }

    private float mOpenScore;
    private float mExtendScore;
}
