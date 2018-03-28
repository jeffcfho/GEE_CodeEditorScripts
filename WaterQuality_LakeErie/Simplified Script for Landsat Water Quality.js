/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var hdpoly = /* color: 0B4A8B */ee.Geometry.Polygon(
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
/// Simplified Script for Landsat Water Quality
// Produces a map of an algal bloom in Lake Erie on 2011/9/3
// Created on 12/7/2015 by Jeff Ho

// VISUALIZATION PARAMETERS -------------------------------------------------------------------------
var startDate = ee.Date.fromYMD(2011,9,3);
//var startDate = ee.Date.fromYMD(2011,8,18);  // has some clouds
//var startDate = ee.Date.fromYMD(2010,9,23);  // very cloudy
var num_days = 16;  // the number of days in the date range

var truecolor = 1;  // show true color image as well
var testThresh = false; // add a binary image classifying into "bloom"and "non-bloom
var bloomThreshold = 0.02346; //threshold for classification fit based on other data
var greenessThreshold = 1.6;

/// FUNCTIONS --------------------------------------------------------------------------------------

// Implements the Automatic Cloud Cover Assessment, with some changes
// http://landsathandbook.gsfc.nasa.gov/pdfs/ACCA_SPIE_paper.pdf
var calcACCA = function(image) {
  var f1 = image.select('B3').lt(0.08); // B3 below 0.08 is non-cloud
  var f2 = image.normalizedDifference(['B2','B5']).gt(0.7); // (B2-B5)/(B2+B5) above 0.7 is non-cloud
  var f3 = image.select('B6').gt(300); // B6 above this is non-cloud
  var f4 = image.expression("( 1 - b('B5') ) * b('B6')").gt(225); // (1-B5)*B6 above this is ambiguous
  var f5 = image.expression("b('B4') / b('B3')").gt(2.0); // B4/B3 above this is ambiguous
  var f6 = image.expression("b('B4') / b('B2')").gt(2.0); // B4/B2 above this is ambiguous
  var f7 = image.expression("b('B4') / b('B5')").lt(1.0); // B4/B5 below this is ambiguous 
  
  // Note: snow not expected in these summer scenes, so filter 8 not used
  //var f8 = image.expression("b('b5') / b('b6')").gt(210); // B5/B6 above this is warm cloud (below is cold cloud)
  
  //ACCA Pass One
  var nonclouds = f1.or(f2).or(f3);
  var ambiguous = nonclouds.not().and(f4.or(f5).or(f6).or(f7));
  var clouds = nonclouds.not().and(ambiguous.not());
  
  return image.addBands(clouds.rename('cloud'));
};

var maskClouds = function(image) {
  // Get regions with clouds:
  var yesCloud = calcACCA(image).select("cloud").eq(0);
  // Remove clouds by making them 0:
  return image.updateMask(yesCloud);
};

// Specifies a threshold for hue to estimate "green" pixels
// this is used as an additional filter to refine the algorithm above
var calcGreenness = function (img) {
  // map r, g, and b for more readable algebra below
  var r = img.select(['B3']);
  var g = img.select(['B2']);
  var b = img.select(['B1']);
  
  // calculate intensity, hue
  var I = r.add(g).add(b).rename(['I']);
  var mins = r.min(g).min(b).rename(['mins']);
  var H = mins.where(mins.eq(r),
    (b.subtract(r)).divide(I.subtract(r.multiply(3))).add(1) );
  H = H.where(mins.eq(g),
    (r.subtract(g)).divide(I.subtract(g.multiply(3))).add(2) );
  H = H.where(mins.eq(b),
    (g.subtract(b)).divide(I.subtract(b.multiply(3))) );
    
  //pixels with hue below 1.6 more likely to be bloom and not suspended sediment
  var Hthresh = H.lte(1.6);
  
  return H.rename('H');
};

// Apply bloom detection algorithm
var calcAlgorithm1 = function(image) {
  // Algorithm 1 based on:
  // Wang, M., & Shi, W. (2007). The NIR-SWIR combined atmospheric 
  //  correction approach for MODIS ocean color data processing. 
  //  Optics Express, 15(24), 15722–15733.
  
  // Add secondary filter using greenness function below
  image = image.addBands(calcGreenness(image));

  // Apply algorithm 1: B4 - 1.03*B5
  var bloom1 = image.select('B4').subtract(image.select('B5').multiply(1.03)).rename('bloom1');

  // Get binary image by applying the threshold
  var bloom1_mask = image.select("H").lte(greenessThreshold).rename(["bloom1_mask"]);
  
  //return original image + bloom, bloom_thresh
  return image.addBands(image)
          .addBands(bloom1)
          .addBands(bloom1_mask);
};

// Apply bloom detection algorithm
var calcAlgorithm2 = function(image) {
  // Algorithm 2 based on: 
  // Matthews, M. (2011) A current review of empirical procedures 
  //  of remote sensing in inland and near-coastal transitional 
  //  waters, International Journal of Remote Sensing, 32:21, 
  //  6855-6899, DOI: 10.1080/01431161.2010.512947
  
  // Apply algorithm 2: B2/B1
  var bloom2 = image.select('B2')
               .divide(image.select('B1'))
               .rename(['bloom2']);
  
  //return original image + bloom, bloom_thresh
  return image.addBands(image)
          .addBands(bloom2);
};


/// APPLY ALGORITHM TO IMAGE -----------------------------------------------------------------------

// Filter to specific image
var collection = ee.ImageCollection('LT5_L1T_TOA')
  .filterDate(startDate, startDate.advance(num_days, 'day'));

// Apply algorithms
collection = collection.map(calcAlgorithm1);
collection = collection.map(calcAlgorithm2);

// Mask out clouds
collection = collection.map(maskClouds);


/// VISUALIZATION --------------------------------------------------------------------

// Color palettes for bloom and bloom_thresh bands
var bPal = ["7C7062", //land -2
            "FFFFFF", //clouds -1
            "0000ff", //non-bloom 0
            "00FF00"  //bloom 1
          ];
var cPal = ['000000', // black 0.014 and below
            'ff00ff', // purple 0.017
            '0000ff', // blue 0.020
            '00ffff', // cyan 0.023
            '00ff00', // green 0.026
            'ffff00', // yellow 0.029
            'ffa500', // orange 0.032
            'ff0000', // red 0.035
            '660000', // red 0.038
            ];

// Create land mask using MODIS annual land cover
var land_mask = ee.Image('MCD12Q1/MCD12Q1_005_2001_01_01').select('Land_Cover_Type_1').neq(0);
var water_mask = ee.Image('MCD12Q1/MCD12Q1_005_2001_01_01').select('Land_Cover_Type_1').eq(0);

// Add layers to map for visualization
Map.addLayer(collection, {bands:['B3','B2','B1'], min:0, max:0.3, 'gamma': 1.2},'True Color');
Map.addLayer(land_mask.mask(land_mask), {min:0, max:1, palette: "000000"}, 'Land Mask');
Map.addLayer(water_mask.mask(water_mask), {min:0, max:1, palette: "AAAAAA"}, 'Water Mask', false);

var bloom_image = collection.mosaic();
Map.addLayer(
  bloom_image.select("bloom1").updateMask(water_mask),
  {min:0.014,max:0.038,palette: cPal},
  'Algorithm 1: Raw Bloom',
  false);
Map.addLayer(
  bloom_image.select("bloom1").updateMask(water_mask.and(bloom_image.select('bloom1_mask'))),
  {min:0.014,max:0.038,palette: cPal},
  'Algorithm 1: Bloom Plus Filter');
Map.addLayer(
  bloom_image.select("bloom2").updateMask(water_mask),
  {min:0.76,max:1,palette: cPal},
  'Algorithm 2: Raw Bloom',
  true);
