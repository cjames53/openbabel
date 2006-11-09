/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.30
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBEdgeBase extends OBBase {
  private long swigCPtr;

  protected OBEdgeBase(long cPtr, boolean cMemoryOwn) {
    super(net.sourceforge.openbabelJNI.SWIGOBEdgeBaseUpcast(cPtr), cMemoryOwn);
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBEdgeBase obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      net.sourceforge.openbabelJNI.delete_OBEdgeBase(swigCPtr);
    }
    swigCPtr = 0;
    super.delete();
  }

  public void setVisit(boolean value) {
    net.sourceforge.openbabelJNI.OBEdgeBase_Visit_set(swigCPtr, this, value);
  }

  public boolean getVisit() {
    return net.sourceforge.openbabelJNI.OBEdgeBase_Visit_get(swigCPtr, this);
  }

  public OBEdgeBase() {
    this(net.sourceforge.openbabelJNI.new_OBEdgeBase__SWIG_0(), true);
  }

  public OBEdgeBase(OBNodeBase arg0, OBNodeBase arg1) {
    this(net.sourceforge.openbabelJNI.new_OBEdgeBase__SWIG_1(OBNodeBase.getCPtr(arg0), arg0, OBNodeBase.getCPtr(arg1), arg1), true);
  }

  public OBGraphBase GetParent() {
    long cPtr = net.sourceforge.openbabelJNI.OBEdgeBase_GetParent(swigCPtr, this);
    return (cPtr == 0) ? null : new OBGraphBase(cPtr, false);
  }

  public void SetParent(OBGraphBase arg0) {
    net.sourceforge.openbabelJNI.OBEdgeBase_SetParent(swigCPtr, this, OBGraphBase.getCPtr(arg0), arg0);
  }

  public long GetIdx() {
    return net.sourceforge.openbabelJNI.OBEdgeBase_GetIdx(swigCPtr, this);
  }

  public void SetIdx(int idx) {
    net.sourceforge.openbabelJNI.OBEdgeBase_SetIdx(swigCPtr, this, idx);
  }

  public void SetBgn(OBNodeBase n) {
    net.sourceforge.openbabelJNI.OBEdgeBase_SetBgn(swigCPtr, this, OBNodeBase.getCPtr(n), n);
  }

  public void SetEnd(OBNodeBase n) {
    net.sourceforge.openbabelJNI.OBEdgeBase_SetEnd(swigCPtr, this, OBNodeBase.getCPtr(n), n);
  }

  public void SwapEnds() {
    net.sourceforge.openbabelJNI.OBEdgeBase_SwapEnds(swigCPtr, this);
  }

  public OBNodeBase GetBgn() {
    long cPtr = net.sourceforge.openbabelJNI.OBEdgeBase_GetBgn(swigCPtr, this);
    return (cPtr == 0) ? null : new OBNodeBase(cPtr, false);
  }

  public OBNodeBase GetEnd() {
    long cPtr = net.sourceforge.openbabelJNI.OBEdgeBase_GetEnd(swigCPtr, this);
    return (cPtr == 0) ? null : new OBNodeBase(cPtr, false);
  }

  public void Error(int f) {
    net.sourceforge.openbabelJNI.OBEdgeBase_Error(swigCPtr, this, f);
  }

  public void SetClosure() {
    net.sourceforge.openbabelJNI.OBEdgeBase_SetClosure(swigCPtr, this);
  }

  public boolean IsAromatic() {
    return net.sourceforge.openbabelJNI.OBEdgeBase_IsAromatic(swigCPtr, this);
  }

  public boolean IsInRing() {
    return net.sourceforge.openbabelJNI.OBEdgeBase_IsInRing(swigCPtr, this);
  }

  public boolean IsClosure() {
    return net.sourceforge.openbabelJNI.OBEdgeBase_IsClosure(swigCPtr, this);
  }

  public boolean Eval(OBEdgeBase arg0) {
    return net.sourceforge.openbabelJNI.OBEdgeBase_Eval(swigCPtr, this, OBEdgeBase.getCPtr(arg0), arg0);
  }

  public long GetBO() {
    return net.sourceforge.openbabelJNI.OBEdgeBase_GetBO(swigCPtr, this);
  }

}