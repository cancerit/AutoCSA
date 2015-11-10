package uk.ac.sanger.cgp.autocsa.util;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;
import org.apache.commons.io.IOUtils;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import uk.ac.sanger.cgp.autocsa.exceptions.CheckedRuntimeCSAException;

/**
 * A class which encapsulates common acitivties against property files, system
 * properties and resource bundles
 * 
 * @author Original: kr2
 * @author  $Author$
 * @version $Revision$
 */
public class PropertyUtils {
  
  protected static final Log log = LogFactory.getLog(PropertyUtils.class.getName());

  /**
   * For a given location this method will load the file into a properties
   * object
   * 
   * @param classpathLocation The loation of the properties file on the 
   * classpath as defined by {@link Class#getResourceAsStream(String)}
   * @return The loaded properties bundle
   * @throws CsaRuntimeException Thrown if an IOException is detected 
   * whilst loading the properties file
   */
  public static Properties getClasspathProperties(String classpathLocation) {
    Properties output = null;
    InputStream is = null;
    try {
      
      if(log.isDebugEnabled()) log.debug("Attempting to load "+classpathLocation);
      
      Class loaderClass = PropertyUtils.class;
      is = new BufferedInputStream(loaderClass.getResourceAsStream(classpathLocation));
      output = new Properties();
      output.load(is);
      
      if(log.isDebugEnabled()) log.debug("Loaded library " + classpathLocation);
    }
    catch(IOException e) {
      throw new CheckedRuntimeCSAException("Could not load library " + classpathLocation, e);
    }
    finally {
      if(is != null) {
        IOUtils.closeQuietly(is);
      }
    }
    return output;
  }
  
  public static Properties getProperties(String configFileLocation) {
    Properties newProperties = null;
    File requestedInput = new File(configFileLocation);
    if(requestedInput.exists()) {
      newProperties = new Properties();
      InputStream is = null;
      try {
        is = new FileInputStream(configFileLocation);
        newProperties.load(is);
      }
      catch(FileNotFoundException e) {
        throw new CheckedRuntimeCSAException("Configuration file could not be found at "+ requestedInput.toString(), e);
      }
      catch(IOException e) {
        throw new CheckedRuntimeCSAException("Configuration file could not be read at "+ requestedInput.toString(), e);
      }
      finally {
        if(is != null) {
          IOUtils.closeQuietly(is);
        }
      }
    }
		return newProperties;
	}
}
