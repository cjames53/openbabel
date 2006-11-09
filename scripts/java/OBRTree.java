/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.30
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBRTree {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected OBRTree(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBRTree obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      net.sourceforge.openbabelJNI.delete_OBRTree(swigCPtr);
    }
    swigCPtr = 0;
  }

  public OBRTree(OBAtom arg0, OBRTree arg1) {
    this(net.sourceforge.openbabelJNI.new_OBRTree(OBAtom.getCPtr(arg0), arg0, OBRTree.getCPtr(arg1), arg1), true);
  }

  public int GetAtomIdx() {
    return net.sourceforge.openbabelJNI.OBRTree_GetAtomIdx(swigCPtr, this);
  }

  public void PathToRoot(SWIGTYPE_p_std__vectorTOpenBabel__OBNodeBase_p_t arg0) {
    net.sourceforge.openbabelJNI.OBRTree_PathToRoot(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBNodeBase_p_t.getCPtr(arg0));
  }

}