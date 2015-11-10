package uk.ac.sanger.cgp.autocsa.exceptions;
                                                                                   
/**
  * New expceiton for using with unchecked exceptions like
  * ArrayIndexOutOfBounds.
  *
  * @author Andy Yates
  * @author $Author$
  * @version $Revision$
  */
public class CheckedRuntimeCSAException extends RuntimeException {
                                                                                   
        public CheckedRuntimeCSAException() {
                super();
        }
                                                                                   
        public CheckedRuntimeCSAException(String msg) {
                super(msg);
        }
                                                                                   
        public CheckedRuntimeCSAException(Throwable t) {
                super(t);
        }
                                                                                   
        public CheckedRuntimeCSAException(String msg, Throwable t) {
                super(msg,t);
        }
}
