package uk.ac.sanger.cgp.autocsa.exceptions;

public class BadTraceException extends Exception {

       public BadTraceException() {
              super();
       }

       public BadTraceException(String message) {
              super(message);
       }

       public BadTraceException(Throwable t) {
              super(t);
       }

       public BadTraceException(String message, Throwable t) {
              super(message, t);
       }

}
