package uk.ac.sanger.cgp.autocsa.beans;

import java.io.Serializable;

public class Mutation implements java.io.Serializable, Comparable {

  private long id_mutation; // can be long rather than Long because primary keys will never be null.
  private int mut_type;
  private int sts_mut_start;
  private int sts_mut_stop;
  private int mut_status;
  private int zygosity;
  private String mut_allele_seq;
  private String wt_allele_seq;
  private boolean flaggedGood = false;
  private String annotationString = "";

  public Mutation(int start, int stop, int type, int status, int zyg, String mut_seq, String wt_seq) {
    sts_mut_start=start;
    sts_mut_stop=stop;
    mut_type=type;
    mut_status=status;
    zygosity=zyg;
    mut_allele_seq=mut_seq;
    wt_allele_seq=wt_seq;
  }
  public Mutation() {
  }
    

    public void setIdMutation(long id)
  {
	  id_mutation = id;
  }
  
  public long getIdMutation()
  {
	  return id_mutation;
  } 
  
  public void setStsMutationStart(int start) {
    sts_mut_start=start;
  }
  public int getStsMutationStart() {
    return sts_mut_start;
  }
  public void setStsMutationStop(int stop) {
    sts_mut_stop=stop;
  }
  public int getStsMutationStop() {
    return sts_mut_stop;
  }
  public void setMutationType(int type) {
    mut_type=type;
  }
  public int getMutationType() {
    return mut_type;
  }
  public void setMutationStatus(int status) {
    mut_status=status;
  }
  public int getMutationStatus() {
    return mut_status;
  }
  public void setZygosity(int zyg) {
    zygosity=zyg;
  }
  public int getZygosity() {
    return zygosity;
  }
  public void setMutationSeq(String seq) {
    mut_allele_seq=seq;
  }
  public String getMutationSeq() {
    return mut_allele_seq;
  }
  public void setWildTypeSeq(String seq) {
    wt_allele_seq=seq;
  }
  public String getWildTypeSeq() {
    return wt_allele_seq;
  }
  
  public boolean isFlaggedGood() {
    return flaggedGood;
  }

  public void setFlaggedGood(boolean flaggedGood) {
    this.flaggedGood = flaggedGood;
  }
	
	public String toString() {
     String nl = System.getProperty("line.separator");
     StringBuffer sb = new StringBuffer();

     String dec = "========";

     sb.append(dec).append(" [S] uk.ac.sanger.cgp.autocsa.beans.Mutation ").append(dec).append(nl);
     sb.append("id_mutation=").append(id_mutation).append(nl);
     sb.append("mut_type=").append(mut_type).append(nl);
     sb.append("sts_mut_start=").append(sts_mut_start).append(nl);
     sb.append("sts_mut_stop=").append(sts_mut_stop).append(nl);
     sb.append("mut_status=").append(mut_status).append(nl);
     sb.append("zygosity=").append(zygosity).append(nl);
     sb.append("mut_allele_seq=").append(mut_allele_seq).append(nl);
     sb.append("wt_allele_seq=").append(wt_allele_seq).append(nl);

     sb.append(dec).append(" [E] uk.ac.sanger.cgp.autocsa.beans.Mutation ").append(dec);

     return sb.toString();
 }
  
  public int compareTo(Object o) {
		Mutation mut = (Mutation) o;
    int result = 0;
		if(this.getStsMutationStart() > mut.getStsMutationStart()) {
      result = 1;
    }
    else {
      result = -1;
    }
    return result;
	}

  public String getAnnotationString() {
    return annotationString;
  }

  public void setAnnotationString(String annotationString) {
    this.annotationString = annotationString;
  }

}
