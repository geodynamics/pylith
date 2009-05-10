class Memory(object):
  sizeInt    = 4
  sizeDouble = 8
  import distutils.sysconfig
  pointerSize = distutils.sysconfig.get_config_var('SIZEOF_VOID_P')
  if pointerSize == 4:
    sizeArrow = 40 # 32 bit
  elif pointerSize == 8:
    sizeArrow = 56 # 64 bit
  elif pointerSize is None:
    sizeArrow = 0 # Use 0 if can't get estimate of pointer size.
  else:
    raise RuntimeError('Could not determine the size of a pointer')

