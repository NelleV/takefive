AFILES = ay2013/rings_10000_raw.bed ay2013/rings_10000_raw.matrix \
	 ay2013/trophozoites_10000_raw.bed \
	 ay2013/trophozoites_10000_raw.matrix \
	 ay2013/schizonts_10000_raw.bed \
	 ay2013/schizonts_10000_raw.matrix \


all: $(AFILES)

ay2013/%:
	mkdir -p ay2013
	wget http://cbio.ensmp.fr/~nvaroquaux/data/$@ -O $@
