/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var hdpoly = /* color: #0b4a8b */ee.Geometry.Polygon(
        [[[-83.57025146484375, 41.69451375429399],
          [-83.3258056640625, 41.66784681529111],
          [-83.089599609375, 41.591888125847646],
          [-82.9632568359375, 41.505556308050686],
          [-82.85064697265625, 41.51172668948318],
          [-82.54302978515625, 41.99732721531279],
          [-82.606201171875, 42.038137999741075],
          [-82.6995849609375, 42.04221763769271],
          [-82.82867431640625, 42.00549146736489],
          [-82.935791015625, 41.98507887314324],
          [-83.03466796875, 42.03609808254534],
          [-83.1060791015625, 42.05445497999412],
          [-83.1060791015625, 42.14208693993462],
          [-83.07861328125, 42.209257234286305],
          [-83.11019897460938, 42.26768689318153],
          [-83.15963745117188, 42.236175276963394],
          [-83.199462890625, 42.09570676103531],
          [-83.2708740234375, 41.9421911163313],
          [-83.33404541015625, 41.94014812213756],
          [-83.3697509765625, 41.8829177113481],
          [-83.47137451171875, 41.792880552910205],
          [-83.51875305175781, 41.70292500835719]]]);
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// Comparing algal bloom detection algorithms (Ho et al., 2017)

// Implements the Automatic Cloud Cover Assessment, with some changes
// http://landsathandbook.gsfc.nasa.gov/pdfs/ACCA_SPIE_paper.pdf
var calcConsACCA = function(img) {
  var f1 = img.select('B3').lt(0.08); // B3 below 0.08 is non-cloud
  var f2 = img.normalizedDifference(['B2','B5']).gt(0.7); // (B2-B5)/(B2+B5) above 0.7 is non-cloud
  var f3 = img.select('B6').gt(300); // B6 above this is non-cloud
  
  /* Do not implement Pass Two here for simplicity, hence these
  estimates are conservative (i.e., all ambiguous  is cloud) */

  var clouds = (f1.or(f2).or(f3)).not();
  
  return img.addBands(clouds.rename('cloud'));
};

//Specifies a threshold for hue to estimate green pixels
var calcGreenness = function (img) {
  var r = img.select(['B3']);
  var g = img.select(['B2']);
  var b = img.select(['B1']);
  var I = r.add(g).add(b).rename(['I']);
  var mins = r.min(g).min(b).rename(['mins']);
  
  var H = mins.where(mins.eq(r),
    (b.subtract(r)).divide(I.subtract(r.multiply(3))).add(1) );
  H = H.where(mins.eq(g),
    (r.subtract(g)).divide(I.subtract(g.multiply(3))).add(2) );
  H = H.where(mins.eq(b),
    (g.subtract(b)).divide(I.subtract(b.multiply(3))) );
  var Hthresh = H.lte(1.6); //threshold of 1.6 fit as described in Ho et al. (2017)
  
  return Hthresh;
};

// Implements the TOA algorithms after correcting for clouds
var calcAlgorithms = function(img) {
  var img2 = calcConsACCA(img);
  var yesCloud = img2.select("cloud");

  //add algorithm outputs as bands
  var img3 = img2.expression("b('B3')/b('B1')").select(["B3"],["RedToBlue"]);
  img3 = img3.addBands(img2.expression("b('B3')/b('B4')").select(['B3'],['RedToNIR']) );
  img3 = img3.addBands(img2.expression("b('B2')/b('B1')").select(['B2'],['GreenToBlue']) );
  img3 = img3.addBands(img2.expression("( b('B1')-b('B3') )/b('B2')")
                           .select(['B1'],['BlueMinusRedOverGreen']) );
  img3 = img3.addBands(img2.expression("47.7-9.21*(2.9594+1.6203*b('B3raw'))/(2.1572+0.9198*b('B1raw'))+" +
                                       "29.7*(3.5046+0.8950*b('B4raw'))/(2.1572+0.9198*b('B1raw'))-" + 
                                       "118*(3.5046+0.8950*b('B4raw'))/(2.9594+1.6203*b('B3raw'))-" + 
                                       "6.81*(3.1591+1.0111*b('B5raw'))/(2.9594+1.6203*b('B3raw'))+" +
                                       "41.9*(2.8122+1.3984*b('B7raw'))/(2.9594+1.6203*b('B3raw'))-" + 
                                       "14.7*(2.8122+1.3984*b('B7raw'))/(3.5046+0.8950*b('B4raw'))")
                           .select(['constant'],['PhycocyaninDetection']) );
  img3 = img3.addBands(img2.select('B4').select(['B4'],['NIR']) );
  img3 = img3.addBands(img2.expression("b('B4')-b('B5')").select(['B4'],['NIRwithSAC']) );
  
  var img_impnirwithsac = img2.expression("b('B4')-1.03*b('B5')").select(["B4"],["ImprovedNIRwithSAC"]); 
  var gness=calcGreenness(img);
  img_impnirwithsac = img_impnirwithsac.where(gness.eq(0),0);
  img3 = img3.addBands(img_impnirwithsac);

  img3=img3.addBands(img2.expression("b('B4') - b('B3')").select(["B4"],["NIRminusRed"]));
  img3=img3.addBands(img2.expression("(b('B4')-b('B5'))/(b('B3')-b('B5'))")
                         .select(["B4"],["NIRoverRedwithSAC"]));
  img3=img3.addBands(img2.expression("( b('B4')-" +
                                        "(b('B4')-b('B1'))+(b('B1')-b('B5'))*(850-490)/(1650-490)" +
                                      ") / " + 
                                      "( b('B3')-" +
                                        "(b('B3')-b('B1'))+(b('B1')-b('B5'))*(660-490)/(1650-490)" + 
                                      ")").select(["B4"],["NIRoverRedwithBAC"]) );
  img3=img3.addBands(img2.expression("(b('B4')-b('B3'))+0.5*(b('B3')-b('B5'))")
                         .select(["B4"],["CurvatureAroundRed"]));
    
  var thresholds = ee.Image([0.462, 0.841, 5.93, 0.0327, 0.0277, 0.0235, 0.495]); 
  var img4=img3.select(['RedToBlue','GreenToBlue','PhycocyaninDetection',
                        'NIR','NIRwithSAC','ImprovedNIRwithSAC','NIRoverRedwithSAC'])
  img4=img4.gte(thresholds);
  
  // mask out clouds
  img4=img4.mask(yesCloud.not());
  
  return (img.addBands(yesCloud).addBands(img4));
};

// Function takes the MERIS_DN band and transforms into CI
var transformMERIS = function(img) {
  img = img.rename('MERIS_DN');
  var justCI = img.select('MERIS_DN');
  var land = justCI.eq(252);
  var clouds = justCI.eq(253);
  
  justCI = justCI.mask(justCI.neq(252)); //mask land
  
  /* convert to CI based on 
      CI = 10.^((double(DN)-10-1)/(250/2.5)-4)
      (itself based on:
        DN =1+(250/2.5)*(4+LOG10(CI))+10
        see Stumpf et al., 2012
      ) 
  */
  justCI = ee.Image(10).pow(
    justCI.double().subtract(10).subtract(1)
    .divide(250/2.5)
    .subtract(4)
    );
    
  //See Stumpf et al. 2012 for thresh of 0.001 ~= 10^5 cells/mL
  justCI=justCI.gte(0.001); 
  justCI = justCI.where(clouds,-0.01/6);
  
  return (img.addBands(justCI.rename('MERIS_CI'))
             .addBands(land.rename('landmask')));
};

var L5_TOA = ee.ImageCollection('LANDSAT/LT5_L1T_TOA'); 
var L5_DN = ee.ImageCollection('LANDSAT/LT5_L1T')
    .select(["B1","B2","B3","B4","B5","B6","B7"],
      ["B1raw","B2raw","B3raw","B4raw","B5raw","B6raw","B7raw"]); 
var L5_TOA_DN = L5_TOA.combine(L5_DN,true); // Phycocyanin algorithm uses DN, so combine with TOA

var erieWBcenter = ee.Geometry.Point(-83.0,41.8);
var img = L5_TOA_DN
  .filterBounds(erieWBcenter)
  .filterDate(ee.Date.fromYMD(2011,9,3).getRange('day'))
  .first();
var img = ee.Image(img); //manual cast so that EE recognizes it as an image
var algs = calcAlgorithms(img);
// Envisat-1 MERIS data provided by European Space Agency 
//  (http://stanford.edu/~jeffho/images/esa_logo.gif)
var merisDN = ee.Image("users/jeffreyh/MERIS_CI/envisat201124609031543CL3GL1e80CI");
var merisCI = transformMERIS(merisDN);

// Default visualization algorithm.
var algorithmToView = 'ImprovedNIRwithSAC'
// one of ['RedToBlue','GreenToBlue','PhycocyaninDetection',
//         'NIR','NIRwithSAC','ImprovedNIRwithSAC','NIRoverRedwithSAC']
var name_list = [
  'ImprovedNIRwithSAC','NIRwithSAC','GreenToBlue',
  'NIRoverRedwithSAC','RedToBlue','NIR','PhycocyaninDetection',
];

var land = ee.Image('ESA/GLOBCOVER_L4_200901_200912_V2_3').select('landcover').neq(210);
land = land.clip(img.geometry()); // clip global landcover map to image
algs = algs.mask(land.not()); //mask water
var land_visParams = {min:1, max:1, palette:['3A4A50']};
var bloom_visParams = {min:0,max:1, palette:['000000','31a354']};
var agree_visParams = {min: 0, max: 2, palette:['000000','ffffff','ff0000',]};
var maps = [];

// Function to switch algorithms
var selectAlg = function(alg_name) {
  //Update left map
  var alg_img = algs.select(alg_name);
  var layer0 = ui.Map.Layer(alg_img,bloom_visParams,'L5 '+alg_name);
  maps[0].layers().set(0, ui.Map.Layer(land,land_visParams,'Land'));
  maps[0].layers().set(1, layer0);
  
  //Update middle map
  var alg_agree = merisCIclipped.neq(algs.select(alg_name));
  var alg_comm = algs.select(alg_name).gt(merisCIclipped);
  alg_agree = alg_agree.where(alg_comm,2);
  var layer1 = ui.Map.Layer(alg_agree,agree_visParams,'AgreementOmissionCommission');
  maps[1].layers().set(0, ui.Map.Layer(land,land_visParams,'Land'));
  maps[1].layers().set(1, layer1);
  
  // Update panel with histogram
  panel.clear();
  var alg_hist = ui.Chart.image.histogram(alg_agree, hdpoly, 90)
    .setSeriesNames([alg_name])
    .setOptions(options);
  panel.add(panelText);
  panel.add(alg_hist);
  panel.add(checkbox);

};

// Map Left - Algorithm to compare
var mapAselect = ui.Select({
  items: name_list,
  placeholder: 'ImprovedNIRwithSAC',
  style:{width:'200px',textAlign:'center'},
  onChange: selectAlg,
});
var mapA = ui.Map();
mapA.add(mapAselect);
mapA.addLayer(land,land_visParams,'Land');
mapA.addLayer(algs.select(algorithmToView),bloom_visParams,'L5 Improved NIRwithSAC');
mapA.addLayer(hdpoly,{color:'FFFFFF'},'Western basin boundaries',false,0.5);
mapA.setControlVisibility(false);
print(mapA.layers())
maps.push(mapA);

// Map Middle - Agreement vs omission errors vs commission errors
var merisCIclipped = merisCI.select('MERIS_CI').clip(img.geometry());
var agree = merisCIclipped.neq(algs.select(algorithmToView));
var commission = algs.select(algorithmToView).gt(merisCIclipped);
agree = agree.where(commission,2);
var mapDiff = ui.Map();
mapDiff.add(ui.Label(
  'Agreement (black) vs. Omission (white) vs. Commission (red)',
  {width:'250px',textAlign:'center'}
));
mapDiff.addLayer(land,land_visParams,'Land');
mapDiff.addLayer(agree,agree_visParams,'AgreementOmissionCommission');
mapDiff.addLayer(hdpoly,{color:'FFFFFF'},'Western basin boundaries',false,0.5);
//red for commission, white for omission, black for agreement
mapDiff.setControlVisibility(false);
print(mapDiff.layers())
maps.push(mapDiff);

// Map Right - MERIS
var mapB = ui.Map();
mapB.add(ui.Label(
  'MERIS CI reference image',
  {width:'200px',textAlign:'center'}
));
mapB.addLayer(merisCI.select('landmask'),land_visParams,'MERIS CI Land');
mapB.addLayer(merisCI.select('MERIS_CI'),bloom_visParams,'Meris CI lake');
mapB.addLayer(hdpoly,{color:'FFFFFF'},'Western basin boundaries',false,0.5);
mapB.setControlVisibility(false);
maps.push(mapB);

var linker = ui.Map.Linker(maps);

// Enable zooming on the top-left map.
maps[0].setControlVisibility({zoomControl: true});

// Show the scale (e.g. '500m') on the bottom-right map.
maps[2].setControlVisibility({scaleControl: true});

// Create a title.
var title = ui.Label('Comparing algal bloom detection algorithms (Ho et al., 2017)', {
  stretch: 'horizontal',
  textAlign: 'center',
  fontWeight: 'bold',
  fontSize: '24px'
});

// Add a left panel to show comparison statistics
var panel = ui.Panel({style: {width: '250px'}});
// Add a label.
var panelText = ui.Label('Comparison statistics vs. MERIS CI over the western basin:');
panel.add(panelText);
// Add a checkbox showing summary area
var checkbox = ui.Checkbox('Show western basin boundaries', false);
checkbox.onChange(function(checked) {
  maps[0].layers().get(2).setShown(checked);
  maps[1].layers().get(2).setShown(checked);
  maps[2].layers().get(2).setShown(checked);
});
// Add a bar chart for Map diff
var options = {
  fontSize: 10,
  hAxis: {
    // title:'Agreement/Omission/Commission',
    ticks: [{v: 0, f: 'Agreement'},
            {v: 1, f: 'Omission'},
            {v: 2, f: 'Commission'}]
  },
  vAxis: {title: 'Number of pixels'},
  colors: ['31a354'],
  // colors: ['#000000','#0000ff','#ff0000'],
};
// Make the histogram, set the options.
var histogram = ui.Chart.image.histogram(agree, hdpoly, 90)
    .setSeriesNames([algorithmToView])
    .setOptions(options);
// I tried to do this to have the histogram colors match the middle image,
//  but doing so seemed to mess up the ticks labels, so I left the bars
//  blue.
// var agree_vis = agree.mask(agree.eq(0));
// agree_vis = agree_vis.addBands(agree.mask(agree.eq(1)));
// agree_vis = agree_vis.addBands(agree.mask(agree.eq(2)));
// agree_vis = agree_vis.rename(['Agreement','Omission','Comission'])
// var histogram = ui.Chart.image.histogram(agree_vis, hdpoly, 90)
//     .setSeriesNames(['Selected Alg vs. MERIS CI'])
//     .setOptions(options);
    
panel.add(histogram);
panel.add(checkbox);

// Create a grid of maps.
var mapGrid = ui.Panel([
    panel,
    ui.Panel([maps[0]], null, {stretch: 'both'}),
    ui.Panel([maps[1]], null, {stretch: 'both'}),
    ui.Panel([maps[2]], null, {stretch: 'both'}),
  ],
  ui.Panel.Layout.Flow('horizontal'), {stretch: 'both'}
);

// Add the maps and title to the ui.root.
ui.root.widgets().reset([title, mapGrid]);
ui.root.setLayout(ui.Panel.Layout.Flow('vertical'));

// Center the maps at Lake Erie.
maps[0].centerObject(erieWBcenter,8);
