/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.30
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class vectorRing {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected vectorRing(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(vectorRing obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      net.sourceforge.openbabelJNI.delete_vectorRing(swigCPtr);
    }
    swigCPtr = 0;
  }

  public vectorRing() {
    this(net.sourceforge.openbabelJNI.new_vectorRing__SWIG_0(), true);
  }

  public vectorRing(long n) {
    this(net.sourceforge.openbabelJNI.new_vectorRing__SWIG_1(n), true);
  }

  public long size() {
    return net.sourceforge.openbabelJNI.vectorRing_size(swigCPtr, this);
  }

  public long capacity() {
    return net.sourceforge.openbabelJNI.vectorRing_capacity(swigCPtr, this);
  }

  public void reserve(long n) {
    net.sourceforge.openbabelJNI.vectorRing_reserve(swigCPtr, this, n);
  }

  public boolean isEmpty() {
    return net.sourceforge.openbabelJNI.vectorRing_isEmpty(swigCPtr, this);
  }

  public void clear() {
    net.sourceforge.openbabelJNI.vectorRing_clear(swigCPtr, this);
  }

  public void add(OBRing x) {
    net.sourceforge.openbabelJNI.vectorRing_add(swigCPtr, this, OBRing.getCPtr(x), x);
  }

  public OBRing get(int i) {
    return new OBRing(net.sourceforge.openbabelJNI.vectorRing_get(swigCPtr, this, i), false);
  }

  public void set(int i, OBRing x) {
    net.sourceforge.openbabelJNI.vectorRing_set(swigCPtr, this, i, OBRing.getCPtr(x), x);
  }

}