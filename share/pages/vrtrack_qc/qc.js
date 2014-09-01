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

// function to sort an observable array of nodes by a property
var sortNodesByProperty = function(list, property, selectlist, selectid) {
    // (we go through this convoluted mess instead of just
    // directly sorting the incoming list, because otherwise the
    // selects don't update until the NEXT call to this method!)
    var arr = [];
    for (var i = list().length - 1; i >= 0; i--) {
        arr.push(list()[i]);
    };
    arr.sort(function(left, right) {
        return left.properties[property] == right.properties[property] ? (left.id < right.id ? -1 : 1) : (left.properties[property] < right.properties[property] ? -1 : 1)
    });
    var selected = [];
    if (selectlist != undefined) {
        for (var i = 0; i < selectlist().length; i++) {
            selected.push(selectlist()[i]);
        };
    }
    list.removeAll();
    list(arr);
    if (selected.length) {
        $(selectid).multiselect('select', selected);
        selectlist(selected);
    }
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
        
        if (method == 'labels') {
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
            getQCGraphData('nodes_of_label', { label: 'Group' }, { resultStore: subargs['groupNodes'], sortProperty: 'name' }, loading, errors);
            fillProperties(labelProperties, 'Study', subargs['studyProperties']);
            fillProperties(labelProperties, 'Donor', subargs['donorProperties']);
            fillProperties(labelProperties, 'Sample', subargs['sampleProperties']);
        }
        else if (method == 'nodes_of_label') {
            resultStore.removeAll();
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
            
            if (subargs.hasOwnProperty('sortProperty') && data.length > 1) {
                sortNodesByProperty(resultStore, subargs['sortProperty']);
            }
        }
        else if (method == 'node_by_id') {
            var result = subargs['result'];
            result(data);
        }
        else if (method == 'donor_qc') {
            var donorGenderResults = subargs['donorGenderResults'];
            var donorInternalDiscordance = subargs['donorInternalDiscordance'];
            donorGenderResults.removeAll();
            donorInternalDiscordance.removeAll();
            var genderArr = [];
            var discArr = [];
            for (var i = 0; i < data.length; i++) {
                var result = data[i];
                var type = result['type'];
                if (type == 'gender') {
                    genderArr.push(result);
                }
                else if (type == 'discordance') {
                    discArr.push(result);
                }
            }
            donorGenderResults(genderArr);
            donorInternalDiscordance(discArr);
        }
    }, subargs);
}