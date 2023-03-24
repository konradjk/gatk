/*
 * $Id: Worklist.java 49076 2007-10-02 20:40:41Z tsharpe $
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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * This class defines a framework for multithreaded processing of a worklist. To use it:<br>
 * 1) Create a class that implements Worklist.<br>
 * 2) Create a class that implements Processor.<br>
 * 3) Create a main method that instantiates an Executive, and tells it to process, using some specified number of
 * threads, an instance of your Worklist subclass using an instance of your Processor subclass.
 * <p>
 * The Executive will spawn the requested number of threads, cloning a Processor for each, and will divvy out items from
 * the Worklist to each Processor, until the Worklist is exhausted.
 * 
 * @author TSharpe
 * @version $Revision: 48897 $
 */
public interface Worklist<T>
{
    /**
     * Get the next work item (represent by a returned object of type T).
     * Return null when the work is done and you want the worker threads to exit.
     * Throwing a WorklistException adds the exception to a list of exceptions that can be
     * inspected when the executive finishes processing the worklist.
     * 
     * @return The next work item (or null to signal the end of the list).
     * @throws WorklistException When it can't get the next work item for some reason.
     */
    T nextWorkitem() throws WorklistException;

    public static class WorklistException
        extends ExceptionWrapper
    {
        /**
         * 
         */
        public WorklistException( String msg, Object workitem )
        {
            super(msg);
            mWorkitem = workitem;
        }

        public WorklistException( String msg, Object workitem, Exception exception )
        {
            super(msg, exception);
            mWorkitem = workitem;
        }

        public Object getWorkitem()
        {
            return mWorkitem;
        }

        private Object mWorkitem;
        private static final long serialVersionUID = 7120952520431866339L;
    }

    /**
     * Makes a worklist from an iterator.
     *
     * @author tsharpe
     * @version $Revision: 48897 $
     */
    public static class IteratorAdaptor<T1>
        implements Worklist<T1>
    {
        public IteratorAdaptor( Iterator<T1> itr )
        {
            if ( itr == null )
            {
                throw new IllegalArgumentException("iterator cannot be null.");
            }
            mItr = itr;
        }

        public synchronized T1 nextWorkitem()
        {
            return mItr.hasNext() ? mItr.next() : null;
        }

        private Iterator<T1> mItr;
    }

    /**
     * Makes a worklist from an array.
     *
     * @author tsharpe
     * @version $Revision: 48897 $
     */
    public static class ArrayAdaptor<T1>
        implements Worklist<T1>
    {
        public ArrayAdaptor( T1[] arr )
        {
            if ( arr == null )
            {
                throw new IllegalArgumentException("array cannot be null.");
            }
            mArr = arr;
            mMaxIdx = arr.length;
        }

        public synchronized T1 nextWorkitem()
        {
            return mIdx < mMaxIdx ? mArr[mIdx++] : null;
        }

        private T1[] mArr;
        private int mMaxIdx;
        private int mIdx;
    }

    /**
     * Makes a worklist by reading a stream line by line.
     *
     * @author tsharpe
     * @version $Revision: 48897 $
     */
    public static class StreamAdaptor
        implements Worklist<String>
    {
        public StreamAdaptor( InputStream strm )
        {
            if ( strm != null )
            {
                mRdr = new BufferedReader(new InputStreamReader(strm));
            }
        }

        public StreamAdaptor( Reader rdr )
        {
            if ( rdr != null )
            {
                mRdr = new BufferedReader(rdr);
            }
        }

        public synchronized String nextWorkitem()
            throws WorklistException
        {
            try
            {
                String result = null;
                if ( mRdr != null )
                {
                    result = mRdr.readLine();
                    if ( result == null )
                    {
                        mRdr.close();
                        mRdr = null;
                    }
                }
                return result;
            }
            catch ( IOException ioe )
            {
                throw new WorklistException("Can't read next workitem from input stream.", null, ioe);
            }
        }

        private BufferedReader mRdr;
    }

    /**
     * Something that processes a work item. Processors must be clonable (so that the Executive can clone one for each
     * thread). The process method does whatever is necessary to handle a work item.
     *
     * @author TSharpe
     * @version $Revision: 48897 $
     */
    interface Processor<T1>
        extends Cloneable
    {
        public Processor<T1> clone();

        /**
         * Process a work item. You'll know what kind of Object to expect as a work item, because you wrote the Worklist
         * subclass to operate in tandem with this one.
         *
         * @param obj
         *            A work item.
         * @throws WorklistException
         *             When the work item can't be processed for some reason.
         */
        public void process( T1 obj ) throws WorklistException;

        /**
         * A signal that processing is complete.
         * @param completedNormally True if the worklist has been processed to exhaustion.
         */
        public void done( boolean completedNormally );
    }

    /**
     * A class to exercise a Processor by retrieving items from the Worklist, one by one, and handing them to the
     * Processor to process. This class is a type of Thread, which exits when the Worklist is exhausted, or when it
     * throws an exception.  (If the processor throws an exception, it's just kept on a list, and the loop continues.)
     *
     * @author TSharpe
     * @version $Revision: 48897 $
     */
    public static class Runner<T1>
        extends Thread
    {
        public Runner( Worklist<T1> worklist, Processor<T1> processor, List<WorklistException> exceptions )
        {
            mWorklist = worklist;
            mProcessor = processor;
            mExceptions = exceptions;
        }

        @Override
        public void run()
        {
            T1 workitem = null;
            boolean complete = false;

            try
            {
                while ( (workitem = mWorklist.nextWorkitem()) != null )
                {
                    try
                    {
                        mProcessor.process(workitem);
                        mNItemsProcessed += 1;
                    }
                    catch ( WorklistException wle )
                    {
                        synchronized ( mExceptions )
                        {
                            mExceptions.add(wle);
                        }
                    }
                    catch ( Exception e )
                    {
                        synchronized ( mExceptions )
                        {
                            mExceptions.add(new WorklistException("Thread " + getName()
                                    + " failed to process workitem '" + workitem + "' due to an unchecked exception.",
                                    workitem, e));
                        }
                    }
                }
                complete = true;
            }
            catch ( Exception e )
            {
                synchronized ( mExceptions )
                {
                    mExceptions.add(new WorklistException("Thread " + getName()
                            + " exited on failure to retrieve the next workitem.", null, e));
                }
            }
            finally
            {
                mProcessor.done(complete);
            }
        }

        public int getNItemsProcessed()
        {
            return mNItemsProcessed;
        }

        Worklist<T1> mWorklist;
        Processor<T1> mProcessor;
        List<WorklistException> mExceptions;
        int mNItemsProcessed;
    }

    /**
     * This is the class you use to kick off the whole can of worms. Its process method spawns a specified number of
     * Runner threads, and waits until all the Runners have exited. When the process method returns, you can gather a
     * variety of entertaining statistics about what happened.
     *
     * @author TSharpe
     * @version $Revision: 48897 $
     */
    public static class Executive
    {
        public Executive()
        {
            mExceptions = new LinkedList<WorklistException>();
        }

        /**
         * Process the worklist in nThreads threads, cloning a processor for each thread.
         *
         * @param nThreads      The number of threads to run.
         * @param worklist      The worklist.
         * @param processor     A processor exemplar.
         */
        @SuppressWarnings("unchecked")
        public <T1> void process( int nThreads, Worklist<T1> worklist, Processor<T1> processor )
        {
            long startMillis = System.currentTimeMillis();

            List<Runner<T1>> threads = new ArrayList<>(nThreads);
            for ( int iii = 0; iii < nThreads; ++iii )
            {
                Runner<T1> thread = new Runner<T1>(worklist, processor.clone(), mExceptions);
                thread.setName("Runner" + iii);
                threads.add(thread);
                thread.start();
            }

            mNItemsProcessed = 0;
            mExceptions.clear();
            for ( final Runner<T1> thread : threads )
            {
                while ( true )
                {
                    try
                    {
                        thread.join();
                        break;
                    }
                    catch ( InterruptedException e )
                    {
                        // nothing to do
                    }
                }

                mNItemsProcessed += thread.getNItemsProcessed();
            }

            long endMillis = System.currentTimeMillis();
            mDuration = endMillis - startMillis;
        }

        public void dumpExceptions()
        {
            Iterator<WorklistException> itr = mExceptions.iterator();
            while ( itr.hasNext() )
            {
                WorklistException wle = itr.next();
                Object workItem = wle.getWorkitem();
                if ( workItem != null )
                    System.err.println("Error: Workitem '" + workItem + "' was not processed.");
                wle.printStackTrace();
            }
        }

        public int getNItemsProcessed()
        {
            return mNItemsProcessed;
        }

        public List<WorklistException> getExceptions()
        {
            return mExceptions;
        }

        public long getDurationMillis()
        {
            return mDuration;
        }

        private int mNItemsProcessed;
        private List<WorklistException> mExceptions;
        private long mDuration;
    }
}
