# XAM2 - experimental XAM rewrite
This repo contains an experimental successor library to XAM.jl.
It is currently work in progress, and not fit for use.

Compared to XAM.jl, it's:
* Far more performance oriented and therefore faster
* Built on a better stack: BGZFLib, CIGARString, XAMAuxData, StringViews, BufferIO
* Has closer adherence to the SAM specs
* Has harder-to-misuse APIs. For example, it does not have uninitialized records,
  is much more explicit about when it can throw, and encodes invariants in types
  to a larger degree.
