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
  
  // For speed of computation, Pass two cloud mask is not used, as would require histogramming each image
  // In this application, precise cloud estimates are not necessary
  // Therefore, to be conservative, ambiguous pixels are treated as clouds
  // Meaning that only the first three thresholds are used
  
  //ACCA Pass One
  var nonclouds = f1.or(f2).or(f3);
  var ambiguous = nonclouds.not().and(f4.or(f5).or(f6).or(f7));
  var clouds = nonclouds.not().and(ambiguous.not());
  
  //ACCA Pass Two
  // No check for desert or snow because not expected in these scenes
  var meanTemp = image.select('B6').mask(clouds).reduceRegion({
        reducer: ee.Reducer.mean(),
        bestEffort: true});
  var p975 = ee.Image(ee.Number(
        image.select('B6').mask(clouds).reduceRegion({
          reducer: ee.Reducer.percentile([97.5]),
          bestEffort: true}).get('B6')
        ));
  var p835 = ee.Image(ee.Number(
        image.select('B6').mask(clouds).reduceRegion({
          reducer: ee.Reducer.percentile([83.5]),
          bestEffort: true}).get('B6')
        ));
  
  var amb_temps = image.select('B6').mask(ambiguous);
  var toAreaKM = ee.Image(30*30/1000/1000);
  var L5_TM_footprint = ee.Image(170*183); //http://landsat.usgs.gov/band_designations_landsat_satellites.php
  
  var uppThermalEffect = amb_temps.lt(p975).and(amb_temps.gt(p835));
  
  var upp_portion = ee.Image(ee.Number(
        uppThermalEffect.reduceRegion({
          reducer: ee.Reducer.sum(),
          bestEffort: true}).get('B6')
        )).multiply(toAreaKM).divide(L5_TM_footprint);
  var upp_meantemp = ee.Image(ee.Number(
        amb_temps.reduceRegion({
          reducer: ee.Reducer.mean(),
          bestEffort: true}).get('B6')
        ));
  
  var lowThermalEffect = amb_temps.lt(p835);
  var low_portion = ee.Image(ee.Number(
        lowThermalEffect.reduceRegion({
          reducer: ee.Reducer.sum(),
          bestEffort: true}).get('B6')
        )).multiply(toAreaKM).divide(L5_TM_footprint);
  var low_meantemp = ee.Image(ee.Number(
        amb_temps.reduceRegion({
          reducer: ee.Reducer.mean(),
          bestEffort: true}).get('B6')
        ));
  
  // Both thermal effect classes become cloud pixels if upper mean < 295K AND
  //  upper thermal effect portion of scene < 0.4
  var ambClouds = (uppThermalEffect.or(lowThermalEffect))
    .multiply(upp_meantemp.lt(295)).multiply(upp_portion.lt(0.4));
  
  // Only lower thermal effect class becomes cloud pixels if low mean < 295K AND
  //  lower thermal effect portion of scene < 0.4
  var ambClouds2 = lowThermalEffect
    .multiply(low_meantemp.lt(295)).multiply(low_portion.lt(0.4));
  
  // Unite ambiguous cloud designations with previous clouds
  clouds = clouds.or(ambClouds.unmask()).or(ambClouds2.unmask());
  
  return image.addBands(clouds.rename('cloud'));
};

//Dates to try: 8/25/2008, 9/10/2008, 9/7/2010, 9/23/2010, 8/9/2011
var TOA_col = ee.ImageCollection('LANDSAT/LT5_L1T_TOA').filterBounds(ee.Feature.Point(-83,41.76))
  //.filterDate(new Date('9/22/2010'),new Date('9/24/2010')); //thin wispy
  .filterDate(new Date('8/24/2008'),new Date('8/26/2008')); //dotty
  //.filterDate(new Date('9/9/2008'),new Date('9/11/2008')); //disrupting scene
var TOA = calcACCA(ee.Image(TOA_col.first()));

Map.addLayer(TOA,{bands:['B3', 'B2', 'B1']}, 'TOA reflectance');
Map.addLayer(TOA.select('cloud'),{},'Cloud');
Map.setCenter(-82.5,41.8,9);