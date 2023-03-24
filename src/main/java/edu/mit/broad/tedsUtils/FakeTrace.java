/**
 * $Id: FakeTrace.java 72836 2008-09-19 15:56:29Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package edu.mit.broad.tedsUtils;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import edu.mit.broad.tedsUtils.Reading.Call;
import edu.mit.broad.tedsUtils.Reading.Sample;

/**
 * Generates a fake trace from sequence.
 * Like sudophred, except works like a filter:  reads sequence from stdin, writes scf to stdout.
 *
 * @author tsharpe
 * @version $Revision: 72836 $
 */
public class FakeTrace
{
    public static void main( String[] args )
        throws IOException
    {
        RefReader rr = new RefReader();
        String seq = rr.getRef();
        createTrace(seq).writeSCF(System.out);
    }

    public static Reading createTrace( String seq )
    {
        int seqLen = seq.length();
        int idealPeakLen = IDEAL_PEAK.length-1;
        Sample[] samples = new Sample[seqLen*idealPeakLen+1];
        ArrayList<Call> calls = new ArrayList<Call>(seqLen);
        int idx = 0;
        double prevA = 0.;
        double prevC = 0.;
        double prevG = 0.;
        double prevT = 0.;
        double weight0 = IDEAL_PEAK[0];
        int callIdx = -idealPeakLen/2;
        for ( int iii = 0; iii < seqLen; ++iii )
        {
            char code = Character.toUpperCase(seq.charAt(iii));
            double fA = 0.;
            double fC = 0.;
            double fG = 0.;
            double fT = 0.;
            switch ( code )
            {
            case 'A': fA = 1.; break;
            case 'C': fC = 1.; break;
            case 'G': fG = 1.; break;
            case 'T': fT = 1.; break;
            }
            calls.add(new Call(code,callIdx+=idealPeakLen,(int)fA*59,(int)fC*59,(int)fG*59,(int)fT*59,0,0,0));
            samples[idx++] = new Sample(weight0*Math.max(fA,prevA),weight0*Math.max(fC,prevC),weight0*Math.max(fG,prevG),weight0*Math.max(fT,prevT));
            prevA = fA;
            prevC = fC;
            prevG = fG;
            prevT = fT;
            for ( int jjj = 1; jjj < idealPeakLen; ++jjj )
            {
                double weight = IDEAL_PEAK[jjj];
                samples[idx++] = new Sample(weight*fA,weight*fC,weight*fG,weight*fT);
            }
        }
        samples[idx] = new Sample(weight0*prevA,weight0*prevC,weight0*prevG,weight0*prevT);
        return new Reading(Arrays.asList(samples),calls,null,null);
    }

    static final double[] IDEAL_PEAK = {32., 128., 320., 464., 576., 624., 640., 624., 576., 464., 320., 128., 32.};
}
