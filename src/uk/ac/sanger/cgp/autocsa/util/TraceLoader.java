package uk.ac.sanger.cgp.autocsa.util;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import org.biojava.bio.chromatogram.Chromatogram;
import org.biojava.bio.chromatogram.ChromatogramFactory;
import org.biojava.bio.chromatogram.UnsupportedChromatogramFormatException;

/**
 * <P>
 * This class is used to load traces into Chromatogram interfaces using the
 * <CODE>ChromatogramFactory.create()</CODE> function. This class should be
 * used if input files are ever going to be gzipped as the class is aware of 
 * this file format and will decompress and create the Chromatogram object
 * automatically.
 * </P>
 *
 * <P>
 * All Chromatogram specific classes can be found in the package <CODE>
 * org.biojava.bio.chromatogram</CODE>
 * </P>
 *
 * <P>
 * Example usage:
 * <BR>
 * <BR>
 * <CODE>
 * File input = new File ("PY834920.scf.gz");<BR>
 * Chromatogram chromat = null;<BR>
 * try {<BR>
 *  chromat = TraceLoader.loadTrace(input);<BR>
 * }<BR>
 * catch (IOException e) {<BR>
 * //NORMAL EXCEPTION CODE<BR>
 * }<BR>
 * catch (UnsupportedChromatogramFormatException e) {<BR>
 * //NORMAL EXCEPTION CODE<BR>
 * }<BR>
 * </CODE><BR>
 * </P><BR>
 *
 * <P>
 * <BR>
 * <CODE>
 * InputStream input = <SCF INPUT STREAM>;<BR>
 * Chromatogram chromat = null;<BR>
 * try {<BR>
 *  chromat = TraceLoader.loadTrace(input);<BR>
 * }<BR>
 * catch (IOException e) {<BR>
 * //NORMAL EXCEPTION CODE<BR>
 * }<BR>
 * catch (UnsupportedChromatogramFormatException e) {<BR>
 * //NORMAL EXCEPTION CODE<BR>
 * }<BR>
 * </CODE><BR>
 * </P>
 *
 * @version $Revision$
 * @author $Author$
 */

public class TraceLoader {
 
  /* 
  This was shown to be 8B1F using od -c | more in unix however this is 
  represented in java as 1F8B for some unknown and bizarre reason I can't be 
  arsed to look up. Fact of the matter is that this works so end of story.
   */
   
  /** Static short representation of the GZip magic number */
  public static final short GZ_MAGIC = (short)0x1F8B;
  
  /**  
   * Takes a trace input, checks to see if it is compressed as .gz (via the 
   * file's magic number) and runs the correct process to load the trace.
   *
   * @param input The input file which is an ABI, SCF or gz compressed versions
   *              of these file types.
   * @return The chromatogram object representation
   * @throws IOException Thrown when an error occurs with the file reading 
   * @throws UnsupportedChromatogramFormatException Thrown when an unknown 
   *            format has been given to the method or that the file is not 
   *              a Trace file.
   */
  public static Chromatogram loadTrace(File input) 
    throws IOException, UnsupportedChromatogramFormatException {
    
    Chromatogram chromat = null;
    
    FileInputStream fis = null;
    
    try {
      fis = new FileInputStream(input);
      chromat = TraceLoader.loadTrace(fis);
    }
    finally {
      fis.close();
    }
    
    return chromat;
    
  }
  
  /**  
   * Takes a trace input, checks to see if it is compressed as .gz (via the 
   * file's magic number) and runs the correct process to load the trace. The
   * method closes the InputStream once finished.
   *
   * @param input The InputStream which is an ABI, SCF or gz compressed versions
   *              of these file types.
   * @return The chromatogram object representation
   * @throws IOException Thrown when an error occurs with the file reading 
   * @throws UnsupportedChromatogramFormatException Thrown when an unknown 
   *            format has been given to the method or that the file is not 
   *              a Trace file.
   */
  public static Chromatogram loadTrace(InputStream input)
    throws IOException, UnsupportedChromatogramFormatException {
    
    Chromatogram chromat = null;
    
    BufferedInputStream bis = new BufferedInputStream(input);
    bis.mark(5);
    DataInputStream dis = new DataInputStream(bis);
    
    //Reading in the first integer in the file which is the magic number
    short out = dis.readShort();
    
    //Nulling dis
    dis = null;

    //Resetting input to pass through to the next section   
    bis.reset();
    
    //Testing to see if this is a GZ file & if so creating the correct input 
    //stream  
    if (out == GZ_MAGIC) {
      BufferedInputStream bufferedSCFGzInput = new BufferedInputStream(new GZIPInputStream(bis)); 
      chromat = ChromatogramFactory.create(bufferedSCFGzInput);
      bufferedSCFGzInput = null;
    }
    //Otherwise this must be an uncompressed trace file
    else {
      chromat = ChromatogramFactory.create(bis);
    }
    
    bis = null;
    
    return chromat;
  }

}