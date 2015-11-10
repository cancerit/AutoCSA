package uk.ac.sanger.cgp.autocsa.exceptions;

public class BadCommentException extends Exception {

       public BadCommentException() {
              super();
       }

       public BadCommentException(String message) {
              super(message);
       }

       public BadCommentException(Throwable t) {
              super(t);
       }

       public BadCommentException(String message, Throwable t) {
              super(message, t);
       }

}
