package uk.ac.sanger.cgp.autocsa.util;

public class AbiIndex {
  public String tag;
  public int list;
  public short dataType;
  public int numItem;
  public long offSet;

  public AbiIndex(String s,int ilist,short dattyp, int num, long oset) {
    tag=s;
    list=ilist;
    dataType=dattyp;
    numItem=num;
    offSet=oset;
  }
}
