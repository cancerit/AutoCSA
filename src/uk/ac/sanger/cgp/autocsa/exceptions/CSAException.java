package uk.ac.sanger.cgp.autocsa.exceptions;

public class CSAException extends Exception {

       public CSAException() {
              super();
       }

       public CSAException(String message) {
              super(message);
       }

       public CSAException(Throwable t) {
              super(t);
       }

       public CSAException(String message, Throwable t) {
              super(message, t);
       }

}
