package uk.ac.sanger.cgp.autocsa.util;

import java.io.File;
import java.io.FilenameFilter;
import java.util.Hashtable;

public class TraceFileFilter implements FilenameFilter {


    private Hashtable filters = null;

    /**
     * Creates a file filter. If no filters are added, then all
     * files are accepted.
     */
    public TraceFileFilter() {
  this.filters = new Hashtable();
    }

    public TraceFileFilter(String well ) {
  this();
  if(well != null) addWell(well);
    }

    public void addWell(String well) {
       if(filters == null) {
           filters = new Hashtable(2);
       }
       filters.put(well, this);
    }

    public String getWell(String filename) {
       if(filename != null) {
         String filenameId=filename.substring(0,2);
         if( filenameId.equals("MP") ) {
           int indexRun=filename.indexOf("Run");
           if( indexRun > 0 ) {
             return filename.substring(indexRun-4,indexRun-1);
           }
         }
         if( filenameId.equals("st") ) {
           int index=filename.length();
           if( index > 0 ) {
//           String ww=filename.substring(index-11,index-8);
//           return filename.substring(index-11,index-8);
           int indexRun=filename.indexOf("Run");
           if( indexRun > 0 ) {
             return filename.substring(indexRun-4,indexRun-1);
           }
           }
         }
       }
       return null;
    }
    public String getWell(File f) {
       if(f != null) {
           String filename = f.getName();
           String result=this.getWell(filename);
           return result;
       }
       return null;
    }

    /**
     * Return true if this file should be shown in the directory pane,
     * otherwise false.
     *
     * Files that begin with "." are ignored.
     */
    public boolean accept(File dir, String name) {
       if(name != null) {
           String well=getWell(name);
           if(well != null && filters.get(well) != null) {
              return true;
           }
       }
       return false;
    }
    public boolean accept(File f) {
  if(f != null) {
      if(f.isDirectory()) {
    return false;
      }
           String well=getWell(f);
           if(well != null && filters.get(well) != null) {
    return true;
      }
  }
  return false;
    }

}
