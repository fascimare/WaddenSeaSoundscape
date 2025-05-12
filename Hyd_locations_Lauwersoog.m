sitepos = [53.41452, 6.23646;... %Lau_W_OFFR
    53.41502, 6.2395;... %Lau_W_ONR
    53.41595, 6.24464;... %Lau_E_OFFR
    53.41615, 6.25273];  %Lau_E_ONR

fykepos = [53.41562, 6.240286;... %SF1
    53.41483,6.238789;... %SF2
    53.4163, 6.248125;... %SF3
    53.41598, 6.246294;... %SF4
    53.4159, 6.251673;... %SF5
    53.4158, 6.253492]; %SF6

%Lat = [53.41452; 53.41502; 53.41595; 53.41615];
%Long = [6.23646; 6.2395; 6.24464; 6.25273];
Lat = [53.41452; 53.41502; 53.41595; 53.41615; 53.41518;  53.41483;  53.4163; 53.41598; 53.4159; 53.4158];
Long = [6.23646; 6.2395; 6.24464; 6.25273; 6.240286; 6.238789; 6.248125; 6.246294; 6.251673; 6.253492];
loc_type = ["Hydrophone";"Hydrophone";"Hydrophone";"Hydrophone"; "SF1"; "SF2"; "SF3"; "SF4"; "SF5"; "SF6"];
Sites = table(Lat,Long,loc_type);
Sites.loc_type = categorical(Sites.loc_type);
lat_lims = [53.414 53.4165];
lon_lims = [6.2 6.3];


h = figure(1);
gb = geobubble(Sites,'Lat','Long', 'BubbleWidthRange',10,'ColorVariable','loc_type','FontSize', 14,'Basemap','satellite');
%gb = geobubble(sitepos(:,1),sitepos(:,2),'ColorVariable',loc_type);
geolimits(lat_lims,lon_lims)
%gb.Basemap = 'landcover';
