RAW DATA structure on aind/scratch

```
📂###### (6 digit mouse ID)/
 ┣ 🖼️vasculature_map.tif
 ┣ 🖼️vasculature_map_annotated.tif* (constantly updated)
 ┣ 📦######_YYYY-MM-DD_HH-MM-SS (for stacks, Session ID by Mouse ID and datetime)/ <- this is the folder that will become a data asset on CodeOcean
 ┃ ┣ 📜rig.json
 ┃ ┣ 📜session.json
 ┃ ┣ 📜subject.json (created during upload to S3, not on VAST)     
 ┃ ┣ 📜data_description.json (created during upload to S3, not on VAST)
 ┃ ┣ 📜metadata.nd.json (created during upload to S3, not on VAST)
 ┃ ┣ 📜procedures.json (created during upload to S3, not on VAST)
 ┃ ┣ 📜processing.json (created during upload to S3, not on VAST)
 ┃ ┃
 ┃ ┗ 📂slap2/
 ┃   ┣ 🖼️vasculature_map_annotated.tif (copied from * upon upload to VAST)
 ┃   ┣ 🖼️session_vasculature_1p.tif
 ┃   ┗ 📂static_data
 ┃     ┣ 📜structure_YYYYMMDD_HHMMSS_DMD#.meta
 ┃     ┣ 💾structure_YYYYMMDD_HHMMSS_DMD#.dat
 ┃     ┣ 🖼️structure_YYYYMMDD_HHMMSS_DMD#.tif
 ┃     ┗ 🖼️structure_YYYYMMDD_HHMMSS_DMD#-REFERENCE.tif
 ┃
 ┣ 📦######_YYYY-MM-DD_HH-MM-SS (for experiments, Session ID by Mouse ID and datetime)/ <- this is the folder that will become a data asset on CodeOcean
 ┃ ┣ 📜rig.json
 ┃ ┣ 📜session.json
 ┃ ┣ 📜subject.json (created during upload to S3, not on VAST)     
 ┃ ┣ 📜data_description.json (created during upload to S3, not on VAST)
 ┃ ┣ 📜metadata.nd.json (created during upload to S3, not on VAST)
 ┃ ┣ 📜procedures.json (created during upload to S3, not on VAST)
 ┃ ┣ 📜processing.json (created during upload to S3, not on VAST)
 ┃ ┃
 ┃ ┣ 📂behavior/
 ┃ ┃ ┣ 📂VCO1_Behavior.harp
 ┃ ┃ ┃ ┣ 💾Behavior_#.bin (# is the register number)
 ┃ ┃ ┃ ┗ 📜device.yml
 ┃ ┃ ┗ 💾<any stim file name>.csv (useful behavior csv files)
 ┃ ┃
 ┃ ┣ 📂behavior-videos/
 ┃ ┃ ┣ 📂BodyCamera
 ┃ ┃ ┃ ┣ 📹video.mp4 (video.avi for now)
 ┃ ┃ ┃ ┗ 📜metadata.csv (metadata.json for now)
 ┃ ┃ ┣ 📂FaceCamera
 ┃ ┃ ┃ ┣ 📹video.mp4 (video.avi for now)
 ┃ ┃ ┃ ┗ 📜metadata.csv (metadata.json for now)
 ┃ ┃ ┗ 📂EyeCamera
 ┃ ┃   ┣ 📹video.mp4 (video.avi for now)
 ┃ ┃   ┗ 📜metadata.csv (metadata.json for now)
 ┃ ┃
 ┃ ┗ 📂slap2/
 ┃   ┣ 🖼️vasculature_map_annotated.tif (copied from * upon upload to VAST)
 ┃   ┣ 🖼️session_vasculature_1p.tif
 ┃   ┗ 📂dynamic_data
 ┃     ┣ 📜acquisition_YYYYMMDD_HHMMSS_DMD#.meta
 ┃     ┣ 💾acquisition_YYYYMMDD_HHMMSS_DMD#-TRIAL######(-CYCLE######).dat
 ┃     ┣ 🖼️acquisition_YYYYMMDD_HHMMSS_DMD#-TRIAL######(-CYCLE######).tif (only sometimes there, depending on SLAP2 mode)
 ┃     ┣ ...
 ┃     ┗ 📂reference_stack
 ┃       ┣ 📜refStack_YYYYMMDD_HHMMSS_DMD#(_CONFIG#).meta
 ┃       ┣ 💾refStack_YYYYMMDD_HHMMSS_DMD#(_CONFIG#).dat
 ┃       ┣ 🖼️refStack_YYYYMMDD_HHMMSS_DMD#(_CONFIG#).tif
 ┃       ┗ 🖼️refStack_YYYYMMDD_HHMMSS_DMD#(_CONFIG#)-REFERENCE.tif
 ┃
 ┣ 📦######_YYYY-MM-DD_HH-MM-SS_slap2_YYYY-MM-DD_HH-MM-SS (first datetime corresponds to the original session data asset, second datetime corresponds to time of starting processing)
 ┃ ┣ 📜rig.json (carried over from raw data)
 ┃ ┣ 📜session.json (carried over from raw data)
 ┃ ┣ 📜subject.json (carried over from raw data)  
 ┃ ┣ 📜data_description.json (carried over from raw data)
 ┃ ┣ 📜metadata.nd.json (carried over from raw data)
 ┃ ┣ 📜procedures.json (carried over from raw data)
 ┃ ┣ 📜processing.json (updated based on processing pipeline that was run)
 ┃ ┣ 📜qc.json (updated based on processing pipeline that was run)
 ┃ ┃
 ┃ ┣ 📜trialTable.mat
 ┃ ┃
 ┃ ┣ 📂motion_correction
 ┃ ┃ ┣ 📜E#T#DMD#_ALIGNMENTDATA.mat
 ┃ ┃ ┣ 🖼️E#T#DMD#_REGISTERED_DOWNSAMPLED-##Hz.tif
 ┃ ┃ ┗ ...
 ┃ ┃
 ┃ ┣ 📂source_extraction
 ┃ ┃ ┣ 📜ExperimentSummary.mat
 ┃ ┃ ┗ 📜ExperimentSummary.h5
 ┃ ┃
 ┃ ┗ 📂qc
 ┃   ┗ ...TBD...
 ┃
 ┗ 📦######_YYYY-MM-DD_HH-MM-SS_<processed modality>_YYYY-MM-DD_HH-MM-SS (first datetime corresponds to the original session data asset, second datetime corresponds to time of starting processing)
   ┣ 📜rig.json (carried over from raw data)
   ┣ 📜session.json (carried over from raw data)
   ┣ 📜subject.json (carried over from raw data)  
   ┣ 📜data_description.json (carried over from raw data)
   ┣ 📜metadata.nd.json (carried over from raw data)
   ┣ 📜procedures.json (carried over from raw data)
   ┣ 📜processing.json (updated based on processing pipeline that was run)
   ┣ 📜qc.json (updated based on processing pipeline that was run)
   ┃
   ┣ 📂<processing step, i.e. annotations>
   ┃ ┗ ...
   ┃
   ┗ 📂qc
```
