%lat_lims = [53 56.0]; %whole Wadden Sea
lat_lims = [53 53.5];
%lat_lims = [52 60]; %North Sea
%lon_lims = [4 10]; %whole Wadden Sea
lon_lims = [4.6 6.85];
%lon_lims = [3 10]; %North Sea
h = figure(1);
geobubble([],[],'Basemap','satellite')
geolimits(lat_lims,lon_lims)

%%
Lat = [53.35; 52.5; 58.5];
Long = [5.5; 4.0; 5.1];
%loc_type = ["Hydrophone";"Hydrophone";"Hydrophone";"Hydrophone"; "SF1"; "SF2"; "SF3"; "SF4"; "SF5"; "SF6"];
Sites = table(Lat,Long);
%Sites.loc_type = categorical(Sites.loc_type);
lat_lims = [52 60];
lon_lims = [6.2 6.3];
h = figure(1);
gb = geobubble(Sites,'Lat','Long', 'BubbleWidthRange',10,'FontSize', 14,'Basemap','satellite','BubbleColorList',[1 0 0.3]);
%gb = geobubble(sitepos(:,1),sitepos(:,2),'ColorVariable',loc_type);
geolimits(lat_lims,lon_lims)
