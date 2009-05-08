class Memory():
  sizeInt    = 4
  sizeDouble = 8
  import distutils.sysconfig
  pointerSize = distutils.sysconfig.get_config_var('SIZEOF_VOID_P')
  if pointerSize == 4:
    sizeArrow = 40 # 32 bit
  elif pointerSize == 8:
    sizeArrow = 56 # 64 bit
  else:
    raise RuntimeError('Could not determine the size of a pointer')

