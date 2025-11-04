mkdir -p topography
cd topography

#Download
curl -L -o 50_120.zip https://dds.cr.usgs.gov/download/eyJpZCI6ODEzODA3MzE3LCJjb250YWN0SWQiOjI3NTc1ODA4fQ==/
curl -L -o 50_150.zip https://dds.cr.usgs.gov/download/eyJpZCI6ODEzODA4MTEwLCJjb250YWN0SWQiOjI3NTc1ODA4fQ==/
curl -L -o 30_150.zip https://dds.cr.usgs.gov/download/eyJpZCI6ODEzNzY2NDk0LCJjb250YWN0SWQiOjI3NTc1ODA4fQ==/
curl -L -o 30_120.zip https://dds.cr.usgs.gov/download/eyJpZCI6ODEzNzQ0ODMwLCJjb250YWN0SWQiOjI3NTc1ODA4fQ==/

unzip 50_120.zip -d tifs_50_120
unzip 50_150.zip -d tifs_50_150
unzip 30_150.zip -d tifs_30_150
unzip 30_120.zip -d tifs_30_120
