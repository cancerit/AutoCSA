package uk.ac.sanger.cgp.autocsa.util;

import java.io.RandomAccessFile;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.io.File;
import java.util.zip.GZIPInputStream;
import java.io.FileInputStream;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
//import newscf.SCF;
import org.biojava.bio.program.scf.SCF;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.chromatogram.Chromatogram;
import org.biojava.bio.chromatogram.ChromatogramFactory;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.chromatogram.UnsupportedChromatogramFormatException;
import org.biojava.bio.symbol.Alignment;
import org.biojava.bio.symbol.IllegalAlphabetException;
import uk.ac.sanger.cgp.autocsa.exceptions.BadCommentException;
import uk.ac.sanger.cgp.autocsa.util.TraceLoader;
import uk.ac.sanger.cgp.autocsa.exceptions.BadTraceException;

/**
 *<p> Class for trace operations, for example reading.</p>
 *
 *Original author:  emd
 *@author $Author$
 *@version $Revision$
 */
public class Trace  {
	
	protected static Log log = LogFactory.getLog(Trace.class.getName());

public String traceFile;
public String traceDir;
public String AbiType;
public int nChannels;
private Channel[] channels; 
private int chanPoints;
private RandomAccessFile fileObj;
private File genericFileObj;
private static int baselineWindow=251;

//New variable for storing a parsed chromatogram
private Chromatogram chromatogram = null;


/**
 * null constructor
 */
public Trace() {
}
/**
 * Allocates a Trace object
 * @param file the trace file File object
 * @param nChans No of Channels in the trace
 */
public Trace(File file, int nChans) {
  this.traceDir=null;
  this.traceFile=null;
  this.genericFileObj=file;
	this.traceFile=file.getPath();
  this.nChannels=nChans;
}
/**
 * Allocates a Trace object
 * @param inputFile the trace file name (full pathname)
 * @param nChans No of Channels in the trace
 */
public Trace(String inputFile, int nChans) {
  this.traceDir=null;
  this.traceFile=inputFile;
  this.fileObj=null;
  this.nChannels=nChans;
  this.genericFileObj=null;
}
/**
 * Allocates a Trace object
 * @param inputFile the trace file name
 * @param inputPath the trace file path
 * @param nChans No of Channels in the trace
 */
public Trace(String inputFile, String inputPath, int nChans) {
  this.traceDir=inputPath;
  this.traceFile=inputPath + "/" + inputFile;
  this.fileObj=null;
  this.nChannels=nChans;
  this.genericFileObj=null;
}

/**
 * Generates a trace object with an already parsed Chromatogram
 * object. Normally used by HoTCSA
 */
public Trace(Chromatogram chromatogram, int nChans) {
	this.traceDir=null;
  this.traceFile=null;
  this.genericFileObj=null;
  this.nChannels=nChans;
  this.chromatogram = chromatogram;
}

/**
 * Generates a trace object with an already parsed Chromatogram
 * object. Normally used by HoTCSA
 */
public Trace(Chromatogram chromatogram, String file, int nChans) {
	this.traceDir=null;
  this.traceFile=file;
  this.genericFileObj=null;
  this.nChannels=nChans;
  this.chromatogram = chromatogram;
}

/**
 * returns the Channel index
 * @param base the base type (one of A,C,G,T)
 * @return the index of the specified Channel type
 */
public int getChannelIndex(String base) {
  String order=Constants.INPUT_BASE_ORDERING;   // order of channels
  return order.indexOf(base);
}
/**
 * returns the maximum trace intensity over an index range
 * @param index1 the start index of the range
 * @param index2 the end index of the range
 * @return the max trace intensity
 */
  public int getMaxTrace(int index1, int index2){
    int max=-30000;
    for (int j=0; j< nChannels; ++j ) {
      max= Math.max(max,channels[j].getMaxIntensity(index1,index2));
    }
    return max;
  }
/**
 * returns the Channel
 * @param index the index of the required Channel
 * @return the specified Channel
 */
  public Channel getChannel(int index) {
    return channels[index];
  }

/**
 * <p>returns the maximum trace intensity over an index range.</p>
 * This is for Heteroduplex traces only and does not include the size
standard channel.
 * @param index1 the start index of the range
 * @param index2 the end index of the range
 * @return the max trace intensity
 */
public int getGlobalMaxHD(int index1,int index2) {
  int maxI=-30000;

// n.b don't include size std channel
  for( int i=0; i < nChannels-1; ++i ) {
    maxI=Math.max(maxI,channels[i].getMaxIntensity(index1,index2));
  }
  return maxI;
}

/**
 * <p>loads a trace file into memory.</p>
 * This method should be called with inputType 1.
 *
 * If the object was constructed with a known chromatogram it will use
 * that chromatogram to derive this object's internal values from
 *
 * @param inputType the input trace type code
 * @throws UnsupportedChromatogramFormatException specific exception
 * @throws IOException if any other miseclleaneous exceptions
 */
public void loadSCFFile (int inputType) throws IOException, UnsupportedChromatogramFormatException, BadTraceException, BadCommentException {

	Chromatogram scfObj=null;

	if(useChromatogram()) {
		scfObj = getChromatogram();
	}
	else {
		File scfFile;
		// open trace file (if we don't already have file pointer)
		if( this.genericFileObj != null ) {
		  scfFile=this.genericFileObj;
		} else {
		  scfFile= new File(this.traceFile);
		}
		
		GZIPInputStream scfGzInput = null;
		try {
		  if( inputType == 1 ) {
		    scfObj = TraceLoader.loadTrace(scfFile);  // use generic trace loader
		  } else if( inputType == 2 ) {
		    scfObj = SCF.create(scfFile);
		  } else {
		    scfGzInput = new GZIPInputStream( new FileInputStream(scfFile));
		    long offset = 0;
		    scfObj = SCF.create(scfGzInput, offset);
		  }
		} catch ( UnsupportedChromatogramFormatException e ) {
		  if(log.isWarnEnabled()) log.warn("Caught UnsupportedChromatogramFormatException: "+e.getMessage());
		  throw new UnsupportedChromatogramFormatException("UnsupportedChromatogramFormatException from TraceLoader.loadTrace()");
		} catch ( IOException e ) {
		  if(log.isWarnEnabled()) log.warn("Caught IOException: "+e.getMessage());
		  throw new IOException("IOException from TraceLoader.loadTrace()"+e.getMessage());
		} finally {
      if(scfGzInput != null) {
        IOUtils.closeQuietly(scfGzInput);
        scfGzInput = null;
      }
		}
	}
	
	int bits=scfObj.getSignificantBits();
	//if(log.isInfoEnabled()) log.info("No of Bits: "+bits);
	chanPoints=scfObj.getTraceLength();
	//if(log.isInfoEnabled()) log.info("No of Points: "+chanPoints);
	try {
	  int m1=scfObj.getMax(DNATools.g());
	  int m2=scfObj.getMax(DNATools.a());
	  int m3=scfObj.getMax(DNATools.t());
	  int m4=scfObj.getMax(DNATools.c());
	//if(log.isInfoEnabled()) log.info("Max I: "+m1+" "+m2+" "+m3+" "+m4);
	} catch ( IllegalSymbolException e ) {
	  if(log.isWarnEnabled()) log.warn("Caught IllegalSymbolException: "+e.getMessage());
	  throw new IOException("IllegalSymbolException N.B Should NOT happen");
	}
  AbiType="3730";   // assumed!!
  channels = new Channel[this.nChannels];
  for( int i=0; i < this.nChannels; ++i ) {
    if( channels[i] == null ) {
      channels[i] = new Channel() ;
    }
  }
  channels[0].name="G";
  channels[1].name="A";
  channels[2].name="T";
  channels[3].name="C";
	try {
	  channels[0].setDataPoints(scfObj.getTrace(DNATools.g()));
	  channels[1].setDataPoints(scfObj.getTrace(DNATools.a()));
	  channels[2].setDataPoints(scfObj.getTrace(DNATools.t()));
	  channels[3].setDataPoints(scfObj.getTrace(DNATools.c()));
	} catch ( IllegalSymbolException e ) {
	  if(log.isWarnEnabled()) log.warn("Caught IllegalSymbolException: "+e.getMessage());
	  throw new IOException("IllegalSymbolException N.B Should NOT happen");
	}
  if(!useChromatogram()) {
    if( isFileAB1Type() ) chromatogram = convertChromatogram(scfObj.getBaseCalls());
  }

} // End method
/** method to convert (Mobilities, Baseline etc) the curent Trace Chromatogram
 * @return Returns a ABIFChromatogramExtender object with updated channels
 */
public ABIFChromatogramExtender convertChromatogram(Alignment baseCall) throws IOException, BadTraceException, BadCommentException {
  ConvertTrace c = new ConvertTrace();
//if( scfObj != null ) c.setChromatogram(scfObj);
  int ok1=c.setChannels(this.traceFile);
  if( ok1 == ConvertTrace.ERROR ) throw new BadTraceException("Failure in Trace Correction (setChannels)");
  if( ok1 == ConvertTrace.IOERROR ) throw new IOException("IO Errror in Trace Correction (setChannels)");
  c.printInfo();
  int ok2=c.processTrace();
  if( ok2 == ConvertTrace.ERROR ) throw new BadTraceException("Failure in Trace Correction (processTrace)");
  if( ok2 == ConvertTrace.IOERROR ) throw new IOException("IO Errror in Trace Correction (processTrace)");
  c.printInfo();
  Channel[] chans=c.getConvertedChannels(Constants.INPUT_BASE_ORDERING);
  for( int j=0; j < this.nChannels; ++j ) {
    channels[j] = chans[j];
  }
  chanPoints=channels[0].getDataPoints().length;
  ABIFChromatogramExtender newChrom=null;
  try {
    newChrom=c.getCorrectedChromatogram(baseCall);
  } catch ( IllegalSymbolException e ) {
    if(log.isWarnEnabled()) log.warn("Caught IllegalSymbolException: "+e.getMessage());
    throw new IOException("IllegalSymbolException N.B Should NOT happen");
  } catch ( IllegalArgumentException e ) {
    if(log.isWarnEnabled()) log.warn("Caught IllegalArgumentException: "+e.getMessage());
    throw new IOException("IllegalArgumentException N.B Should NOT happen");
  } catch ( IllegalAlphabetException e ) {
    if(log.isWarnEnabled()) log.warn("Caught IllegalAlphabetException: "+e.getMessage());
    throw new IOException("IllegalAlphabetException N.B Should NOT happen");
  }
  return newChrom;
}

public boolean isFileAB1Type(String file) throws FileNotFoundException, IOException {
  int magic = -1;
  RandomAccessFile raf = null;
  try {
    raf = new RandomAccessFile(file,"r");
    magic=raf.readInt();
  }
  catch(FileNotFoundException e) {
    throw e;
  }
  catch(IOException e) {
    throw e;
  }
  finally {
    if(raf != null) {
      raf.close();
    }
  }
  
  if ( magic ==  ChromatogramFactory.ABI_MAGIC ) {
    return true;
  } else {
    return false;
  }
}

private boolean isFileAB1Type() throws FileNotFoundException, IOException {
  return isFileAB1Type(this.traceFile);
}

/**
 * performs a baselining of a specified Channel
 * @param traceIndex index of Channel to baseline
 * @param traceLim1 start index of Channel to baseline
 * @param traceLim2 end index of Channel to baseline
 */
public void baseline( int traceIndex, int traceLim1,int traceLim2) {

/**
; Function; calculate a value for the local baseline at each data
;           point and subtracts from data array, uses a window of
;           width wlen at each data point
;           N.B. wlen must be of form (2n+1)
;           If traceLim1 and traceLim2 are 0 then these are calculated and
;           indicate the usefull range of the trace to use for
;           baselining, i.e to save computation
;           (Calulated using the size standard channel)
;
; Note: This function is intended to apply the same processes to
;       the raw data channels as the auto-analysis process.
; The only 2 relevant data operation that constitute auto-analysis
; are 1) baseline correction
;     2) smoothing
; In the analysis module G5_CSCE.gsp (used operationally on the ABI 3100's)
;     1) baseline window width= 251
;     2) smoothing set to OFF
*/

int hlen=(baselineWindow-1)/2;
int np=this.chanPoints;
int[] data=this.channels[traceIndex].getDataPoints();
int i,i1,i2;

int[] bline=new int[np];
int[] bline1=new int[np];
for ( i=0 ; i<np ; ++i ) {
  bline[i]=0;
  bline1[i]=0;
}
if ( traceLim2 == 0 ) { traceLim2=this.chanPoints-1; }
for ( i=traceLim1 ; i<=traceLim2 ; ++i ) {
  i1=Math.max(i-hlen,0);
  i2=Math.min(i+hlen,this.chanPoints-1);
  bline[i]=this.channels[traceIndex].getMinIntensity(i1,i2);
}
// temporarily set trace data channel to the intermediate baseline bline
// so we can apply the method getMaxIntensity to this data (bline)
this.channels[traceIndex].setDataPoints(bline);

for ( i=traceLim1 ; i<=traceLim2 ; ++i ) {
  i1=Math.max(i-hlen,0);
  i2=Math.min(i+hlen,this.chanPoints-1);
  bline1[i]=this.channels[traceIndex].getMaxIntensity(i1,i2);
}

for ( i=0 ; i<np ; ++i ) {
  data[i]-=bline1[i];
// set any -ve signal to zero
  if ( data[i] < 0 ) { data[i]=0; }
}

// reset original trace data to baselined
this.channels[traceIndex].setDataPoints(data);
}
/**
 * derives a usefull lower limit index for baselining.
 * Note Only to be used for Heteroduplex type traces
 * @return the lower limit to use for baselining
 */
public int getbaselineLimit1 () {
  int hlen=(baselineWindow-1)/2;
  int np=this.chanPoints;

// find pos of dye flare in Size Std Channel
  int[] data=this.channels[4].getDataPoints();
  int indexFlare=0;
  int maxval=0;
  for ( int i=0 ; i<np-1 ; ++i ) {
    if ( data[i+1] > maxval ) {
      maxval=data[i+1];
      indexFlare=i+1;
    }
  }
  if ( indexFlare-100 > hlen ) {
    return indexFlare-100;
  } else {
    return hlen;
  }
}
/**
 * derives a usefull upper limit index for baselining.
 * Note Only to be used for Heteroduplex type traces
 * @return the upper limit to use for baselining
 */
public int getbaselineLimit2 () {
  int hlen=(baselineWindow-1)/2;
  int np=this.chanPoints;

// find pos of last Size Std by checking large -ve gradients
  int[] data=this.channels[4].getDataPoints();
  int indexGrad=0;
  int grad;
  for ( int i=0 ; i<np-1 ; ++i ) {
    grad=data[i+1]-data[i];
    if ( grad < -40 ) { indexGrad=i+1; }
  }
  if ( indexGrad > 0 ) {
    if ( indexGrad+200 < np-1-hlen ) {
      return indexGrad+200;
    } else {
      return np-1-hlen;
    }
  } else {
    return np-1-hlen;
  }
}

	public Chromatogram getChromatogram() {
		return chromatogram;
	}

	/**
	 * Checks the internal objects chromatogram and file_position_from_root
	 * for null values and returns a boolean indicating if it is possible to
	 * use the chromatogram object or not
	 *
	 * @return Indicates if you should use the internal Chromatogram object
	 * or you should retrieve the chromatogram from {@link #getFilePositionFromRoot()}
	 */
	public boolean useChromatogram() {
		boolean useChromatogram = false;
		if(getChromatogram() != null) {
			useChromatogram = true;
		}
		return useChromatogram;
	}
}
