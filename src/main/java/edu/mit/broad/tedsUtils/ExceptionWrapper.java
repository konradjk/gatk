/*
 * $Id: ExceptionWrapper.java 49076 2007-10-02 20:40:41Z tsharpe $
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

import java.io.PrintStream;
import java.io.PrintWriter;

/**
 * Class that handles wrapping other exceptions.
 * @author Ted Sharpe
 * @version $Revision: 48897 $
 */
public class ExceptionWrapper
    extends Exception
{
    /**
     * Vanilla constructor.
     */
    public ExceptionWrapper( String msg )
    {
        super(msg);
    }

    /**
     * Nested-exception constructor.
     */
    public ExceptionWrapper( String msg, Exception exception )
    {
        super(msg, exception);
    }

    @Override
    public void printStackTrace()
    {
        if ( getCause() == null )
            super.printStackTrace();
        else
        {
            System.err.println(this + " Caused by:");
            getCause().printStackTrace();
        }
    }

    @Override
    public void printStackTrace( PrintStream strm )
    {
        if ( getCause() == null )
            super.printStackTrace(strm);
        else
        {
            strm.println(this + " Caused by:");
            getCause().printStackTrace(strm);
        }
    }

    @Override
    public void printStackTrace( PrintWriter wrtr )
    {
        if ( getCause() == null )
            super.printStackTrace(wrtr);
        else
        {
            wrtr.println(this + " Caused by:");
            getCause().printStackTrace(wrtr);
        }
    }

    private static final long serialVersionUID = -3728322112418595716L;
}
