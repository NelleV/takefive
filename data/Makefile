AFILES = ay2013/rings_10000_raw.bed ay2013/rings_10000_raw.matrix \
	 ay2013/trophozoites_10000_raw.bed \
	 ay2013/trophozoites_10000_raw.matrix \
	 ay2013/schizonts_10000_raw.bed \
	 ay2013/schizonts_10000_raw.matrix \

LFILES = lemieux2013/25kb/NF54_combined_raw.bed \
	 lemieux2013/25kb/NF54_combined_raw.matrix \
	 lemieux2013/25kb/DCJ_Off_combined_raw.bed \
	 lemieux2013/25kb/DCJ_Off_combined_raw.matrix \
	 lemieux2013/25kb/DCJ_On_combined_raw.bed \
	 lemieux2013/25kb/DCJ_On_combined_raw.matrix \
	 lemieux2013/25kb/B15C2_combined_raw.bed \
	 lemieux2013/25kb/B15C2_combined_raw.matrix \
	 lemieux2013/25kb/A4_pilot_combined_raw.bed \
	 lemieux2013/25kb/A4_pilot_combined_raw.matrix \
	 lemieux2013/25kb/A44_combined_raw.bed \
	 lemieux2013/25kb/A44_combined_raw.matrix \
	 lemieux2013/25kb/A4+_combined_raw.bed \
	 lemieux2013/25kb/A4+_combined_raw.matrix \


all: $(AFILES) $(LFILES)

ay2013/%:
	mkdir -p ay2013
	wget http://cbio.ensmp.fr/~nvaroquaux/data/$@ -O $@

lemieux2013/%:
	mkdir -p $(@D)
	wget http://cbio.ensmp.fr/~nvaroquaux/data/$@ -O $@
