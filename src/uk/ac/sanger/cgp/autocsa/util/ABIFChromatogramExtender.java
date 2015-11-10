package uk.ac.sanger.cgp.autocsa.util;

import org.biojava.bio.program.abi.ABIFChromatogram;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Alignment;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;

/** ABIFChromatogramExtender Object to return a modified ABIFChromatogram
 * @author emd
 */
public class ABIFChromatogramExtender extends ABIFChromatogram {


/** Working constructor for ABIFChromatogramExtender
 * @param newChans array[4] of Channel objects holding new channel values
 * @param nSigBits int holding No. of required significant bits to hold data
 */
public ABIFChromatogramExtender(Channel[] newChans, int nSigBits, Alignment baseCall) throws IllegalArgumentException, IllegalSymbolException, IllegalAlphabetException {

  setBits(nSigBits);
  setTrace(DNATools.c(),newChans[0].getDataPoints(),newChans[0].getMaxIntensity());
  setTrace(DNATools.a(),newChans[1].getDataPoints(),newChans[1].getMaxIntensity());
  setTrace(DNATools.g(),newChans[2].getDataPoints(),newChans[2].getMaxIntensity());
  setTrace(DNATools.t(),newChans[3].getDataPoints(),newChans[3].getMaxIntensity());
	setBaseCallAlignment(baseCall);
	
	
}
}
