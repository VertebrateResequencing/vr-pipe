function getQueryStringParameterByName(name) {
    name = name.replace(/[\[]/, "\\[").replace(/[\]]/, "\\]");
    var regex = new RegExp("[\\?&]" + name + "=([^&#]*)"),
        results = regex.exec(location.search);
    return results == null ? "" : decodeURIComponent(results[1].replace(/\+/g, " "));
}

// function to fill an array with the properties of a given label
var fillProperties = function(labelProperties, label, list) {
    var props = labelProperties[label];
    var arr = [];
    for (var i = 0; i < props.length; i++) {
        arr.push(props[i]);
    }
    list(arr);
}

// function to call qc methods to get results from the graph db
var getQCGraphData = function(method, args, subargs, loading, errors) {
    loading.removeAll();
    errors.removeAll();
    
    vrpipeRestMethod('qc', method, args, loading, errors, function(data) {
        var resultStore = subargs['resultStore'];
        // we won't push directly to this because even with
        // ratelimiting it's too slow; we'll push to temp arrays
        // and then set the resultStore to it
        
        switch (method) {
            case 'labels':
                resultStore.removeAll();
                var keys = [];
                for (var key in data) {
                    keys.push(key);
                }
                keys.sort();
                var arr = [];
                var labelProperties = subargs['labelProperties'];
                for (var i = 0; i < keys.length; i++) {
                    var key = keys[i];
                    labelProperties[key] = data[key];
                    if (key != 'Group' && key != 'Study') {
                        arr.push(key);
                    }
                }      
                resultStore(arr);
                var viewLabels = subargs['viewLabels'];
                viewLabels.push('Donor');
                viewLabels.push('Sample');
                
                // populate the groups select and the all properties
                fillProperties(labelProperties, 'Study', subargs['studyProperties']);
                fillProperties(labelProperties, 'Donor', subargs['donorProperties']);
                fillProperties(labelProperties, 'Sample', subargs['sampleProperties']);
                getQCGraphData('nodes_of_label', { label: 'Group' }, { resultStore: subargs['groupNodes'], sortProperty: 'name' }, loading, errors);
                break;
            
            case 'nodes_of_label':
                var arr = [];
                var flatten = subargs['flatten'];
                for (var i = 0; i < data.length; i++) {
                    if (flatten) {
                        data[i]['properties']['node_id'] = data[i]['id'];
                        data[i]['properties']['node_label'] = data[i]['label'];
                        arr.push(data[i]['properties']);
                    }
                    else {
                        arr.push(data[i]);
                    }
                }
                resultStore(arr);
                break;
            
            case 'node_by_id':
                var result = subargs['result'];
                result(data);
                break;
            
            case 'donor_qc':
                var donorGenderResults = subargs['donorGenderResults'];
                var donorFluidigmDiscordance = subargs['donorFluidigmDiscordance'];
                var donorGenotypingDiscordance = subargs['donorGenotypingDiscordance'];
                var donorCopyNumberSummary = subargs['donorCopyNumberSummary'];
                var donorAberrantRegions = subargs['donorAberrantRegions'];
                var donorAberrantPolysomy = subargs['donorAberrantPolysomy'];
                var donorCopyNumberPlot = subargs['donorCopyNumberPlot'];
                var donorLOHCalls = subargs['donorLOHCalls'];
                var donorPluritestSummary = subargs['donorPluritestSummary'];
                var donorPluritestPlots = subargs['donorPluritestPlots'];
                donorGenderResults.removeAll();
                donorFluidigmDiscordance.removeAll();
                donorGenotypingDiscordance.removeAll();
                donorCopyNumberSummary.removeAll();
                donorAberrantRegions.removeAll();
                donorAberrantPolysomy.removeAll();
                donorCopyNumberPlot(undefined);
                donorLOHCalls.removeAll();
                donorPluritestSummary.removeAll();
                donorPluritestPlots.removeAll();
                var genderArr = [];
                var discFArr = [];
                var discGArr = [];
                var cnsArr = [];
                var arArr = [];
                var apArr = [];
                var lohArr = [];
                var psArr = [];
                var ppArr = [];
                
                for (var i = 0; i < data.length; i++) {
                    var result = data[i];
                    var type = result['type'];
                    delete result['type'];
                    switch (type) {
                        case 'gender':
                            genderArr.push(result);
                            break;
                        case 'discordance_fluidigm':
                            discFArr.push(result);
                            break;
                        case 'discordance_genotyping':
                            discGArr.push(result);
                            break;
                        case 'copy_number_summary':
                            cnsArr.push(result);
                            break;
                        case 'aberrant_regions':
                            for (var key in result) {
                                if (key == 'graph' && result.hasOwnProperty(key)) {
                                    result[key] = '/file' + result[key];
                                    if (result[key] == '/filenull') {
                                        result[key] = '-none-';
                                    }
                                }
                            }
                            arArr.push(result);
                            break;
                        case 'aberrant_polysomy':
                            for (var key in result) {
                                if (key == 'graph' && result.hasOwnProperty(key)) {
                                    result[key] = '/file' + result[key];
                                }
                            }
                            apArr.push(result);
                            break;
                        case 'copy_number_plot':
                            donorCopyNumberPlot('/file' + result['plot']);
                            break;
                        case 'loh_calls' :
                            lohArr.push(result);
                            break;
                        case 'pluritest_summary' :
                            psArr.push(result);
                            break;
                        case 'pluritest_plot':
                            result['path'] = '/file' + result['path'];
                            ppArr.push(result);
                            break;
                    }
                }
                
                donorGenderResults(genderArr);
                donorFluidigmDiscordance(discFArr);
                donorGenotypingDiscordance(discGArr);
                donorCopyNumberSummary(cnsArr);
                donorAberrantRegions(arArr);
                donorAberrantPolysomy(apArr);
                donorLOHCalls(lohArr);
                donorPluritestSummary(psArr);
                donorPluritestPlots(ppArr);
                break;
            
            case 'sample_discordance':
                var sdResults = subargs['sdResults'];
                sdResults.removeAll();
                var arr = [];
                for (var i = 0; i < data.length; i++) {
                    arr.push(data[i]);
                }
                sdResults(arr);
                break;
            
            default:
                errors.push('invalid qc method: ' + method);
        }
    }, subargs);
}