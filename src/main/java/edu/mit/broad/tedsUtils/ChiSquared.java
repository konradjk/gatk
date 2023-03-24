/*
 * $Id: ChiSquared.java 51142 2007-11-05 17:19:49Z tsharpe $
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
 * A bunch of static methods to calculate chi-squared, and it's cumulative distribution
 * function. Algorithms adapted from "Numerical Recipes in C/The Art of Scientific
 * Computing" by William H. Press, Brian P. Flannery, Saul A. Teukolsky, and William T.
 * Vettering, Cambridge University Press, (c) 1988.
 * @author Ted Sharpe
 * @version $Revision: 23649 $
 */
public class ChiSquared
{
    public static double chiSq( ContingencyTable ds )
    {
        int nX = ds.getNX();
        double[] xTot = new double[nX];
        int nY = ds.getNY();
        double[] yTot = new double[nY];
        double grandTot = 0.;

        for ( int iX = 0; iX < nX; ++iX )
        {
            for ( int iY = 0; iY < nY; ++iY )
            {
                double vvv = ds.get(iX, iY);
                xTot[iX] += vvv;
                yTot[iY] += vvv;
                grandTot += vvv;
            }
        }

        double result = 0.;

        for ( int iX = 0; iX < nX; ++iX )
        {
            for ( int iY = 0; iY < nY; ++iY )
            {
                double exp = xTot[iX] * yTot[iY] / grandTot;
                double diff = ds.get(iX, iY) - exp;
                result += diff * diff / exp;
            }
        }

        return result;
    }

    /**
     * Calculates the probability that a chi squared value as large, or larger, than the
     * one supplied could occur by chance.
     * @param chiSq The chi squared value.
     * @param nu The number of degrees of freedom.
     * @return Probability of chance occurrence.
     */
    public static double probChiSq( double chiSq, double nu )
    {
        return MathFunctions.incompleteGamma(nu / 2., chiSq / 2.);
    }
}
