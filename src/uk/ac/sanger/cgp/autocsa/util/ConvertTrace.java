package uk.ac.sanger.cgp.autocsa.util;

import java.io.File;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.nio.ByteBuffer;
import org.apache.commons.io.IOUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import org.biojava.bio.program.abi.ABIFChromatogram;
import org.biojava.bio.program.abi.ABIFParser;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.chromatogram.AbstractChromatogram;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.chromatogram.UnsupportedChromatogramFormatException;
import org.biojava.bio.symbol.Alignment;
import org.biojava.bio.symbol.IllegalAlphabetException;

import uk.ac.sanger.cgp.autocsa.beans.MobilityCorData;
import uk.ac.sanger.cgp.autocsa.exceptions.BadCommentException;

/** ConvertTrace Object to perform Trace data Channel conversion,
 * i.e. Baselining, Mobility Correction and signal start & end removal
 * @author emd
 */
public class ConvertTrace extends AbstractChromatogram {
	
protected static Log log = LogFactory.getLog(ConvertTrace.class.getName());

private static final String internalBaseOrder="CAGT";
public static final int OK=0;
public static final int ERROR=-1;
public static final int IOERROR=-2;
private static final int baselineWinLen=251;
private int signalStart=0;
private int signalEnd=0;
private int primerLoc=0;
private int numPoints=0;
private int numChans=4;
private float nominalAvgBaseSpac=12.0f;

private static final String ABI_CTAG1="SVER";
private static final String ABI_CTAG2="SPAC";
private static final String ABI_CTAG3="PDMF";
private static final String SCF_CTAG1="VER1";
private static final String SCF_CTAG2="SPAC";
private static final String SCF_CTAG3="DYEP";
private float udcVersion=0.0f;
private float traceAvgBaseSpac=0.0f;
private String dyepLabel=null;
private String mobFileName=null;
//private Chromatogram abiObj=null;
private ABIFChromatogram abiObj=null;
private MobilityCorData[] mobCor=null;
private Channel[] channels=null;

public ConvertTrace() {}

public void setABIFChromatogram(ABIFChromatogram c) {
  this.abiObj=c;
}
public ABIFChromatogram getABIFChromatogram() {
  return this.abiObj;
}
protected AbstractChromatogram reverseComplementInstance() {
  return null;
}
public int processTrace() throws BadCommentException {

/*
** Remove baseline offset.
** n.b. best to baseline first so no edge effects due to baseline window
*/
    if( subtractBaseline() == ERROR )
    {
      if(log.isInfoEnabled()) log.info("processTrace: error: bad status: subtractBaseline");
      return( ERROR );
    }

/*
** Find start of signal.
*/
    if( findSignalStart() == ERROR )
    {
      if(log.isInfoEnabled()) log.info("processTrace: error: bad status: findStart");
      return( ERROR );
    }

/*
** Trim preprimer sequence.
*/
    trimTraceStart(signalStart);

/*
** Trim off end portion of trace
*/
    trimEndTrace();

/*
** Trim primer peak.
*/
    trimTraceStart(primerLoc - signalStart);

/*
** Perform mobility correction.
*/
    constructMobFileName();
    if( mobilityCorrect() == IOERROR )
    {
      if(log.isInfoEnabled()) log.info("processTrace: error: bad status: mobilityCorrect");
      return( IOERROR );
    }

    return( OK );
}
public Channel[] getConvertedChannels() {
// n.b. internal order of Channels is CAGT
// but required in order GATC (same order as AB1/SCF file)
  Channel[] chans = new Channel[4];
  chans[0]=channels[2];
  chans[1]=channels[1];
  chans[2]=channels[3];
  chans[3]=channels[0];
  return chans;
}
public Channel[] getConvertedChannels(String order) {
// n.b. internal order of Channels is CAGT
  if( order.length() != numChans ) return null;
  order=order.toUpperCase();
  Channel[] chans = new Channel[4];
  chans[0]=channels[internalBaseOrder.indexOf(order.substring(0,1))];
  chans[1]=channels[internalBaseOrder.indexOf(order.substring(1,2))];
  chans[2]=channels[internalBaseOrder.indexOf(order.substring(2,3))];
  chans[3]=channels[internalBaseOrder.indexOf(order.substring(3,4))];
  return chans;
}
public void readCommentsBlockABI(File abiFile) throws IOException {
  ABIFParser abiParse = null;
  FileInputStream fis = null;
  try {
    fis = new FileInputStream(abiFile);
    abiParse=new ABIFParser(fis);
    String s=null;
    ABIFParser.DataAccess a1=abiParse.getDataAccess();
    ABIFParser.TaggedDataRecord r1=abiParse.getDataRecord(ABI_CTAG1,1);
    ABIFParser.TaggedDataRecord r2=abiParse.getDataRecord(ABI_CTAG2,1);
    ABIFParser.TaggedDataRecord r3=abiParse.getDataRecord(ABI_CTAG3,1);
  //f[3] =(float)Double.longBitsToDouble( r2.dataRecord); // doesn't work
    traceAvgBaseSpac =Float.intBitsToFloat( (int) r2.dataRecord);
    byte[] str = new byte[8];
    ByteBuffer.wrap(str).putLong(r1.dataRecord);
    udcVersion = Float.parseFloat(new String(str));

  // n.b. this tag is a p-string (first byte is string len) so ignore
    byte[] str1 = new byte[(int)r3.recordLength-1];
    for( int k=1 ; k < r3.recordLength; ++k) str1[k-1]=r3.offsetData[k];
    dyepLabel=new String(str1);
  }
  finally {
    IOUtils.closeQuietly(fis);
  }

}
public void readCommentsBlockSCF() {

/*
  Properties comBlock=scfObj.getComments();
  udcVersion=Float.parseFloat(comBlock.getProperty(SCF_CTAG1));
  traceAvgBaseSpac=Float.parseFloat(comBlock.getProperty(SCF_CTAG2));
  dyepLabel=comBlock.getProperty(SCF_CTAG3);
  if(log.isInfoEnabled()) log.info("TAGS: "+udcVersion+" "+traceAvgBaseSpac+" "+dyepLabel);
*/

}
private void constructMobFileName() throws BadCommentException {
// get the MOBDIR environmental variable
  String mobdir="."+ File.separator + "resources" + File.separator +"mobCorrFiles" + File.separator;

// map DYEP tag value to real mobfilename
// here we change a '{' to a '=' and skip any '}'
  
  if(dyepLabel.trim().length() == 0) {
    if(log.isWarnEnabled()) log.warn("No dyep comment was found, mobility correction cannot be performed");
    throw new BadCommentException("No dyep comment was found, mobility correction cannot be performed");
  }
  
  dyepLabel=dyepLabel.replaceAll("\\{","=");
  dyepLabel=dyepLabel.replaceAll("\\}","");
  mobFileName=mobdir.concat(dyepLabel);

// test if the derived file for the mobility correction file exists
// if not then use a generic version for current machine & chemistry
  if ( ! (new File(mobFileName)).exists() ) {
    String pop=null,machine=null;
    if( dyepLabel.indexOf("POP4") > -1 ) pop="POP4";
    if( dyepLabel.indexOf("POP5") > -1 ) pop="POP5";
    if( dyepLabel.indexOf("POP6") > -1 ) pop="POP6";
    if( dyepLabel.indexOf("POP7") > -1 ) pop="POP7";
    if( dyepLabel.indexOf("3100") > -1 ) machine="3100";
    if( dyepLabel.indexOf("3130") > -1 ) machine="3130";
    if( dyepLabel.indexOf("3700") > -1 ) machine="3700";
    if( dyepLabel.indexOf("3730") > -1 ) machine="3730";
    mobFileName=mobdir.concat("Generic"+machine+pop+".mob");
  }
  if(log.isInfoEnabled()) log.info("MOB: "+mobFileName);
}
private int trimTraceStart(int len) {
  int i, j;

  /*
  ** return if nothing to do.
  */
  if( len == 0 ) {
    return( OK );
  }

  for( j = 0; j < numChans; ++j ) {
    int k=0;
    int[] trace=channels[j].getDataPoints();
    int[] t=new int[numPoints-len];
    for( i = len; i < numPoints; ++i ) {
      t[k++] = trace[i];
    }
    channels[j].setDataPoints(t);
  }
  numPoints -= len;

  return( OK );
}
private int trimEndTrace() {
  int i, j;
  /*
  ** trim off end of trace
  ** n.b signalEnd is original trace pos - subtract signalStart
  */
  if( signalEnd >= numPoints ) return( OK );

  int len=signalEnd-signalStart;
  for( j = 0; j < numChans; ++j ) {
    int[] trace=channels[j].getDataPoints();
    int t[] = new int[len];
    for( i = 0; i < len; ++i ) {
      t[i] = trace[i];
    }
    channels[j].setDataPoints(t);
  }
  numPoints = len;

  return( OK );

}
private int applyMobilityCorrect(int[] trace, int startOffset, MobilityCorData mobCor, int chan) {
  int i, j, k, m;

  int[] scan = mobCor.getScan();
  int[] shift = mobCor.getShift();
  int numShift=mobCor.getNumShift();
  int[] sdata = new int[numPoints];
  /*
  ** perform channel shifting
  */
  j = shift[0];
  k = 1;
  for( i = 0; i < numPoints; ++i ) {
    m = i - startOffset;
    if( k < numShift && m >= scan[k] ) {
      j = shift[k];
      ++k;
    }

    if( i + j < 0 ) {
      sdata[i] = trace[0];
    } else if( i + j >= numPoints ) {
      sdata[i] = trace[numPoints - 1];
    } else {
      sdata[i] = trace[j+i];
    }
  }
  /*
  ** copy values back to trace
  */
  for( i = 0; i < numPoints; ++i ) {
    trace[i] = sdata[i];
  }
  channels[chan].setDataPoints(trace);

  return( OK );
}

private int readMobilityCorrect() {
  int i, j;
  int numChanInFile=-1;
  /*
  ** Open file.
  */
  File iFile=new File(mobFileName);
  BufferedReader in = null;
  String line=null;
  try {
    in = new BufferedReader(new FileReader(iFile));
    // Read number of channels.
    line=in.readLine();
    numChanInFile=Integer.parseInt(line);
    if( numChanInFile < numChans )
    { 
      if(log.isWarnEnabled()) log.warn("readMobilityCorrect: error: number of channels in mobility file");
      if(log.isWarnEnabled()) log.warn("is less than number in trace data "+numChanInFile);
      in.close();
      return( IOERROR );
    }

  /*
  ** Instantiate mobCor array
  */
    mobCor=new MobilityCorData[numChans];
    for( j = 0; j < numChans; ++j ) {
      mobCor[j] = new MobilityCorData();
    }

  /*
  ** Read line that contain the number of correction
  ** factors for each channel.
  */
  line=in.readLine();
  String[] strs=line.split("\\s+");
  if( strs.length < numChans )
  { 
    if(log.isWarnEnabled()) log.warn("readMobilityCorrect: error: corrupted mobility file");
    in.close();
    return( IOERROR );
  }

  /*
  ** Extract numbers of correction factors from line.
  */
  for( j = 0; j < numChans; ++j ) {
    mobCor[j].setNumShift(Integer.parseInt(strs[j]));
  }

  /*
  ** Read correction factors and save to beans
  */
  for( j = 0; j < numChans; ++j ) {
    mobCor[j].setDoCorrection(1);
    int np=mobCor[j].getNumShift();
    int[] scan = new int[np];
    int[] shift = new int[np];
    for( i = 0; i < np; ++i ) {
      line=in.readLine();
      strs=line.split("\\s+");
      int curChan=Integer.parseInt(strs[0]);
      scan[i]= Integer.parseInt(strs[1]);
      shift[i]= Integer.parseInt(strs[2]);
      if( curChan != j ) {
        if(log.isWarnEnabled()) log.warn("readMobilityCorrect: unexpected channel in mobfile");
        in.close();
        return( IOERROR );
      }
    }
    if( scan[np-1] == 0 && shift[np-1] == 0) mobCor[j].setDoCorrection(0);
    mobCor[j].setScan(scan);
    mobCor[j].setShift(shift);
  }

    in.close();
  } catch ( IOException e ) {
    if(log.isWarnEnabled()) log.warn("Caught IOException: "+e.getMessage());
    return( IOERROR );
  }
  finally {
    if(in != null) {
      IOUtils.closeQuietly(in);
    }
  }

  return( OK );
}

private int mobilityCorrect() {
  int j, k;
  float scale=1.0f;
  
  /*
  ** read mobility correction values.
  */
  if( readMobilityCorrect() == IOERROR )
  {
    if(log.isWarnEnabled()) log.warn("mobilityCorrect: problem reading mobfile "+mobFileName);
    return( IOERROR );
  }
  int startOffset = primerLoc - signalStart;

// make sure av base spacing is set zero for UDC 1.0 traces
  if( udcVersion < 2.0f ) traceAvgBaseSpac=0.0f;

  if( traceAvgBaseSpac > 0.0f )
  {
    scale = traceAvgBaseSpac / nominalAvgBaseSpac;
  }
  if(log.isInfoEnabled()) log.info("Mobility scale factor: "+scale);
  
  for( j = 0; j < numChans; ++j ) {
    if( mobCor[j].getDoCorrection() == 1 ) {
      int[] trace=channels[j].getDataPoints();
      int[] scan = mobCor[j].getScan();
      int[] shift = mobCor[j].getShift();
// n.b original C code did not use "round" so truncates floats
// hence using round as below is more accurate
      for ( k = 0; k < mobCor[j].getNumShift(); k++ ) {
        scan[k] = Math.round( (float) scan[k] * scale );
        shift[k] = Math.round( (float) shift[k] * scale );
      }
      mobCor[j].setScan(scan);
      mobCor[j].setShift(shift);
      applyMobilityCorrect(trace, startOffset, mobCor[j], j );
    }
  }

  return( OK );
}
private int findSignalStart() {
  int j, iStart;
// cut off specified scans from each end of the trace to avoid dye blobs
// and larger spikes at end of trace - to get a reasonable max intensity
  int[] minPeakAmps = new int[numChans];
  int cut1=2500;
  int cut2=1000;
  for( j=0; j < numChans; ++j ) {
    minPeakAmps[j]=channels[j].getMaxIntensity(cut1,numPoints-cut2);
    minPeakAmps[j]=(int) (0.04f*(float) minPeakAmps[j]);
    minPeakAmps[j]=Math.max(250,minPeakAmps[j]);
  }
// do peak search (with shoulders turned off)
  int minpos=10000,maxpos=0;
  int nPeaks=0;
  for( j=0; j < numChans; ++j ) {
    channels[j].findPeaks(minPeakAmps[j],40000,40000,5);
    int[] pos=channels[j].getPeaksPos();
    if( pos == null ) continue;
    nPeaks+=pos.length;
    minpos=Math.min(pos[0],minpos);
    maxpos=Math.max(pos[pos.length-1],maxpos);
  }
// if < 50 peaks return an error status
  if( nPeaks < 50 ) return( ERROR );

// choose a value that's the highest multiple of 50
// and that is < minpos (n.b. same principle as original plan code)
  iStart= (int) Math.floor((float) minpos / 50.0f) * 50;
  signalEnd=Math.min(maxpos + 500,numPoints); // add 500 for safety
  signalStart   = ( iStart - 10 > 0 ) ? iStart - 10 : 0;
  primerLoc     = iStart;

  return( OK );
}
private int subtractBaseline() {
  int i, ii, j, k, m;
  float val;

  if( numPoints < baselineWinLen ) {
    if(log.isInfoEnabled()) log.info("subtractBaseline: error - Not enough data points: "+numPoints);
    return( ERROR );
  }

  float temp1[] = new float[numPoints];
  float temp2[] = new float[numPoints];

  for( j = 0; j < numChans; ++j ) {
    int[] trace=channels[j].getDataPoints();
    for( m = 0; m < 2; m++) {
      /*
      ** apply the max min filter m=0 is min m=1 is max
      */
      for( i = 0; i < numPoints; ++i ) {
        if( m == 0 ) {
          temp1[i] = trace[i];
        } else {
          temp1[i] = temp2[i];
        }
      }
      for( i = 0; i < numPoints; ++i ) {
        ii = i + 1;
        if( ii >= baselineWinLen ) {
          val = temp1[ii-baselineWinLen];
          if( m == 1 ) {
            for( k = ii - baselineWinLen; k < ii; ++k ) {
              if( temp1[k] > val ) {
                val = temp1[k];
              }
            }
          } else {
            for( k = ii - baselineWinLen; k < ii; ++k ) {
              if( temp1[k] < val ) {
                val = temp1[k];
              }
            }
          }
          temp2[ii-baselineWinLen/2-1] = val;
        }
      }
      for( i = 0; i < baselineWinLen / 2; ++i ) {
        temp2[i] = temp2[baselineWinLen/2];
        temp2[numPoints-i-1] = temp2[numPoints-baselineWinLen/2-1];
      }
    }

//  Subtract derived baseline
    for( i = 0; i < numPoints; ++i ) {
      trace[i] -= temp2[i];
    }
    channels[j].setDataPoints(trace);
  }

  return( OK );
}
public int setChannels(String tFile) {
  File abiFile= new File(tFile);
  try {
//  scfObj = TraceLoader.loadTrace(scfFile);
    if( abiObj == null ) {
      FileInputStream fis = null;
      try {
        fis = new FileInputStream(abiFile);
        abiObj=ABIFChromatogram.create(fis);
      }
      finally {
        IOUtils.closeQuietly(fis);
      }
      //abiObj=ABIFChromatogram.create(abiFile);
    }
  } catch ( UnsupportedChromatogramFormatException e ) {
    if(log.isWarnEnabled()) log.warn("Caught UnsupportedChromatogramFormatException: "+e.getMessage());
    return(IOERROR);
  } catch ( IOException e1 ) {
    if(log.isWarnEnabled()) log.warn("Caught IOException: "+e1.getMessage());
    return(IOERROR);
  }
  channels = new Channel[this.numChans];
  for( int i=0; i < this.numChans; ++i ) {
    if( channels[i] == null ) {
      channels[i] = new Channel() ;
    }
  }
// n.b. required order of Channels is CAGT
  channels[0].name="C";
  channels[1].name="A";
  channels[2].name="G";
  channels[3].name="T";
  try {
    channels[0].setDataPoints(abiObj.getTrace(DNATools.c()));
    channels[1].setDataPoints(abiObj.getTrace(DNATools.a()));
    channels[2].setDataPoints(abiObj.getTrace(DNATools.g()));
    channels[3].setDataPoints(abiObj.getTrace(DNATools.t()));
    numPoints=channels[0].getChannelLength();
    readCommentsBlockABI(abiFile);
  } catch ( IllegalSymbolException e ) {
    if(log.isWarnEnabled()) log.warn("Caught IllegalSymbolException: "+e.getMessage());
    return(IOERROR);
  } catch ( IOException e1 ) {
    if(log.isWarnEnabled()) log.warn("Caught IOException: "+e1.getMessage());
    return(IOERROR);
  }
  return(OK);
}
public ABIFChromatogramExtender getCorrectedChromatogram(Alignment baseCall) throws IllegalArgumentException, IllegalSymbolException, IllegalAlphabetException {
  int nSigBits=0;
  if( abiObj != null) nSigBits=abiObj.getSignificantBits();
//ABIFChromatogramExtender abiObj=new ABIFChromatogramExtender(channels,nSigBits);
  return new ABIFChromatogramExtender(channels,nSigBits,baseCall);
}
public void printInfo() {
  if(log.isInfoEnabled()) log.info("TRACE: "+numPoints);
  if(log.isInfoEnabled()) log.info("Signal range: "+signalStart+" "+signalEnd);
}
}
