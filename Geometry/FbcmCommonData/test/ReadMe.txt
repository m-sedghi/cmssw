the "fbcmGeometryXML_cfi.py" should be located in the following path

/Geometry/TrackerCommonData/python/fbcmGeometryXML_cfi.py


pwd --> /Geometry/TrackerCommonData/fbcmtest

## No need for Compilation -->>  cd ../ && scram b -j 16 && cd fbcmtest && cmsRun dumpFBCM_cfg.py <<-
>> cmsRun dumpFBCM_cfg.py
>> cmsShow -c simGeo.fwc --sim-geom-file=cmsSimGeom.root
or
>> cmsShow -c unnamed.fwc --sim-geom-file=cmsSimGeom.root
