/*
 * $Id: FastaBit.java 73459 2008-10-10 16:10:03Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2008 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package edu.mit.broad.tedsUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.LinkedList;
import java.util.List;

/**
 * Utility to read segmented FASTA files.
 *
 * @author tsharpe
 */
public class FastaBit
{
    /**
     * Make one.
     */
    public FastaBit( String name, String sequence )
    {
        mName = name;
        mSequence = sequence;
    }

    public String getName()
    {
        return mName;
    }

    public String getSequence()
    {
        return mSequence;
    }

    /**
     * Reads segmented FASTA files, one section at a time.
     */
    public static class Rdr
    {
        /**
         * Make one.
         */
        public Rdr( Reader rawRdr )
        {
            mRdr = new BufferedReader(rawRdr);
        }

        /**
         * Read the next comment-delimited chunk of FASTA.
         * Returns null if there are no more chunks.
         * The reader will be closed if you read until this method returns null.
         * @throws IOException
         */
        public FastaBit next()
            throws IOException
        {
            FastaBit result = null;

            if ( mRdr != null )
            {
                StringBuilder sb = new StringBuilder();
                String line;

                while ( (line = mRdr.readLine()) != null )
                {
                    if ( line.length() > 0 )
                    {
                        if ( line.charAt(0) == '>' ) 
                        {
                            String lastComment = mNextComment;
                            mNextComment = line.substring(1).trim();
                            if ( sb.length() > 0 )
                            {
                                return new FastaBit(lastComment,sb.toString()); // EARLY RETURN!
                            }
                        }
                        else
                        {
                            sb.append(line);
                        }
                    }
                }
                close();
                mRdr = null;
    
                if ( sb.length() > 0 )
                {
                    result = new FastaBit(mNextComment,sb.toString());
                }
            }

            return result;
        }

        /**
         * Call close if you don't read successfully to the end.
         * (Or just close the reader yourself.)
         * @throws IOException
         */
        public void close()
            throws IOException
        {
            mRdr.close();
        }

        private BufferedReader mRdr;
        private String mNextComment;
    }

    /**
     * Reads all the FASTA, closes the reader.
     */
    public static List<FastaBit> create( Reader rawRdr )
        throws IOException
    {
        List<FastaBit> result = new LinkedList<FastaBit>();
        Rdr rdr = new Rdr(rawRdr);
        FastaBit bit;
        while ( (bit = rdr.next()) != null )
        {
            result.add(bit);
        }
        return result;
    }

    private String mName;
    private String mSequence;
}
