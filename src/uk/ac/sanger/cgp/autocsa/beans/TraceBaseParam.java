package uk.ac.sanger.cgp.autocsa.beans;

import java.io.Serializable;

/**
 *  TraceBaseParam - a bean of four primitives. The first (int) contains the position
 in the sts of the base, the second (float) contains the base quality value,
 the third (int) the scan index of the base peak, and fourth (int) the base intensity
 */

public class TraceBaseParam implements java.io.Serializable
{
    
    private int stsPosition;
    private float traceBaseQuality;
    private int traceBaseScanIndex;
    private int traceBaseIntensity;

     /**
     *  Shortcut constructor
     *  @param position STS position of the base.
     *  @param quality float value.
     *  @param scan int value.
     *  @param intensity int value.
     */
    public TraceBaseParam(int position, float quality, int scan, int intensity)
    {
        stsPosition = position;
        traceBaseQuality = quality;
        traceBaseScanIndex = scan;
        traceBaseIntensity = intensity;
    }
    
    /**
     *  Setter for setting all 3 class fields
     *  @param position STS position of the base.
     *  @param quality int value.
     *  @param scan int value.
     *  @param intensity int value.
     */
    public void setTraceParams(int position, float quality, int scan, int intensity)
    {
        stsPosition = position;
        traceBaseQuality = quality;
        traceBaseScanIndex = scan;
        traceBaseIntensity = intensity;
    }
       
    /**
     *  Setter for setting all 3 class fields
     *  @param position STS position of the base.
     *  @param quality int value.
     *  @param scan int value.
     *  @param intensity int value.
     */
    public void setTraceParams(int position, int quality, int scan, int intensity)
    {
        stsPosition = position;
        traceBaseQuality = quality;
        traceBaseScanIndex = scan;
        traceBaseIntensity = intensity;
    }

    /**
     * Gets the current value of stsPosition
     * @return Current value of stsPosition
     */
    public int getPosition()
    {
        return stsPosition;
    }
    
    /**
     * Gets the current value of traceBaseQuality
     * @return Current value of traceBaseQuality
     */
    public float getBaseQuality()
    {
        return traceBaseQuality;
    }

    /**
     * Gets the current value of traceBaseScanIndex
     * @return Current value of traceBaseScanIndex
     */
    public int getBaseScanIndex()
    {
        return traceBaseScanIndex;
    }

    /**
     * Sets the current value of traceBaseScanIndex
     * @param New value of traceBaseScanIndex
     */
    public void setBaseScanIndex(int scanIndex)
    {
        this.traceBaseScanIndex=scanIndex;
    }

    /**
     * Gets the current value of traceBaseIntensity
     * @return Current value of traceBaseIntensity
     */
    public int getBaseIntensity()
    {
        return traceBaseIntensity;
    }
    
    public String toString() {
        String nl = System.getProperty("line.separator");
        StringBuffer sb = new StringBuffer();

        String dec = "========";

        sb.append(dec).append(" [S] uk.ac.sanger.cgp.autocsa.beans.TraceBaseParam ").append(dec).append(nl);
        sb.append("stsPosition=").append(stsPosition).append(nl);
        sb.append("traceBaseQuality=").append(traceBaseQuality).append(nl);
        sb.append("traceBaseScanIndex=").append(traceBaseScanIndex).append(nl);
        sb.append("traceBaseIntensity=").append(traceBaseIntensity).append(nl);

        sb.append(dec).append(" [E] uk.ac.sanger.cgp.autocsa.beans.TraceBaseParam ").append(dec);

        return sb.toString();
    }
}

