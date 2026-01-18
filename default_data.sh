# Get default static field
mkdir -p "$1/gravity/static"
wget -P "$1/gravity/static" https://icgem.gfz-potsdam.de/getmodel/gfc/17de3c1a5391d7d0397425eb7cfc8367ccbf77e33b8709600550ff8ea7ae9ceb/GOCO05c.gfc
wget -P "$1/gravity/static" https://ftp.tugraz.at/outgoing/ITSG/GRACE/ITSG-Grace2018/static/ITSG-Grace2018s.gfc

# Get default GRACE data

# Get default normals data
mkdir -p "$1/gravity/monthly/normals"
wget -P "$1/gravity/monthly/normals" ftp://ftp.tugraz.at/outgoing/ITSG/GRACE/ITSG-Grace_operational/monthly/normals_SINEX/monthly_n96/ITSG-Grace_operational_n96_2024-12.snx.gz

# Get default AOD1B data file
mkdir -p "$1/aod1b"
wget -P "$1/aod1b" ftp://isdcftp.gfz.de/grace/Level-1B/GFZ/AOD/RL07/2025/AOD1B_2025-01-01_X_07.asc.gz

# Decompress any .gz file in data directory
gzip -d -r "$1"