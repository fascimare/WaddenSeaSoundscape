sitepos = [53.41452, 6.23646;... %Lau_W_OFFR
    53.41502, 6.2395;... %Lau_W_ONR
    53.41595, 6.24464;... %Lau_E_OFFR
    53.41615, 6.25273];  %Lau_E_ONR

Lat = [53.41452; 53.41502; 53.41595; 53.41615];
Long = [6.23646; 6.2395; 6.24464; 6.25273];
loc_type = ["Hydrophone";"Hydrophone";"Hydrophone";"Hydrophone"];
Sites = table(Lat,Long,loc_type);
Sites.loc_type = categorical(Sites.loc_type);
lat_lims = [53.414 53.4165];
lon_lims = [6.2 6.3];


h = figure(1);
gb = geobubble(Sites,'Lat','Long', 'BubbleWidthRange',10,'ColorVariable','loc_type','FontSize', 14,'Basemap','satellite');
%gb = geobubble(sitepos(:,1),sitepos(:,2),'ColorVariable',loc_type);
geolimits(lat_lims,lon_lims)
%gb.Basemap = 'landcover';
