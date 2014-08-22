// we have a custom ko binding to handle our loading feedback
// (from http://blog.greatrexpectations.com/2012/06/17/loading-placeholders-using-knockout-js/)
ko.bindingHandlers.loadingWhen = {
    init: function (element) {
        var $element = $(element),
            currentPosition = $element.css("position")
            $loader = $("<div>").addClass("loader").hide();
        
        //add the loader
        $element.append($loader);
        
        //make sure that we can absolutely position the loader against the original element
        if (currentPosition == "auto" || currentPosition == "static")
            $element.css("position", "relative");

        //center the loader
        $loader.css({
            position: "absolute",
            top: "20px",
            left: "50%",
            "margin-left": -($loader.width() / 2) + "px"
        });
    },
    update: function (element, valueAccessor) {
        var isLoading = ko.utils.unwrapObservable(valueAccessor()),
            $element = $(element),
            $childrenToHide = $element.children(":not(div.loader)"),
            $loader = $element.find("div.loader");

        if (isLoading) {
            $childrenToHide.css("visibility", "hidden").attr("disabled", "disabled");
            $loader.fadeIn("slow");
        }
        else {
            $loader.fadeOut("fast");
            $childrenToHide.css("visibility", "visible").removeAttr("disabled");
        }
    }
};

// we can have sortable table columns
// (from http://stackoverflow.com/a/16964843/675083)
ko.bindingHandlers.sort = {
    init: function(element, valueAccessor, allBindingsAccessor, viewModel, bindingContext) {
        var asc = false;
        element.style.cursor = 'pointer';
        
        element.onclick = function(){
            var value = valueAccessor();
            var prop = value.prop;
            var data = value.arr;
            
            asc = !asc;
            
            data.sort(function(left, right){
                var rec1 = left;
                var rec2 = right;
                
                if(!asc) {
                    rec1 = right;
                    rec2 = left;
                }
                
                var props = prop.split('.');
                for(var i in props){
                    var propName = props[i];
                    var parenIndex = propName.indexOf('()');
                    if(parenIndex > 0){
                        propName = propName.substring(0, parenIndex);
                        rec1 = rec1[propName]();
                        rec2 = rec2[propName]();
                    } else {
                        rec1 = rec1[propName];
                        rec2 = rec2[propName];
                    }
                }
                
                return rec1 == rec2 ? 0 : rec1 < rec2 ? -1 : 1;
            });
        };
    }
};

// ajax function to send/get json queries
var ajax = function(uri, method, loading, error, data) {
    var loadingTimeout;
    var request = {
        url: uri,
        type: method,
        contentType: "application/json",
        accepts: "application/json",
        cache: false,
        dataType: 'json',
        data: JSON.stringify(data),
        beforeSend: function (xhr) {
            // say that we're loading, but only if for more than
            // half a second
            if (loadingTimeout) { clearTimeout(loadingTimeout); }
            loadingTimeout = setTimeout(function() {
                loading.push(true);
            }, 500);
        },
        error: function(jqXHR, textStatus, errorThrown) {
            if (textStatus) {
                var msg = "Ajax " + textStatus + " " + jqXHR.status + " connecting to " + uri;
                if (errorThrown) {
                    msg = msg + ": " + errorThrown;
                }
                error.push(msg);
            }
        },
        complete: function(jqXHR, textStatus) {
            if (loadingTimeout) { clearTimeout(loadingTimeout); }
            loading.pop();
        }
    };
    return $.ajax(request);
}

// function to call a VRPipe Rest method and do something with the
// json response, and its helper function to deal with errors
var vrpipeRestDataErrorParser = function(data, error) {
    if (data == null) {
        error.push('no data');
        return 1;
    }
    if (data.hasOwnProperty('errors')) {
        for (var i = 0; i < data.errors.length; i++) {
            error.push(data.errors[i]);
        }
        return 1;
    }
    return 0;
}
var vrpipeRestMethod = function(category, method, args, loading, error, donefunction, other) {
    var promise = ajax('/rest/' + category + '/' + method, 'POST', loading, error, args);
    
    if (donefunction) {
        promise.done(function(data) {
            if (! vrpipeRestDataErrorParser(data, error)) {
                donefunction(data, other);
            }
        });
    }
    else {
        return promise;
    }
}

var shortenString = function(string) {
    var shortened = string;
    if (string.length > 25) {
        shortened = string.substring(0, 5) + '[...]' + string.slice(-15);
    }
    return shortened;
}

var wbr = function(string) {
    var wbrd = '' + string;
    if (string.length > 21) {
        wbrd = wbrd.replace(/(.{21})/g, "$1<wbr>");
    }
    else {
        wbrd = wbrd.replace(/([_\-,\.\(\)])/g, "$1<wbr>");
    }
    return wbrd;
}

var zeropad = function(num) {
    return (num < 10) ? ("0" + num) : num;
}

// convert epoch seconds to friendly date
var epochToDate = function(epoch) {
    var date = new Date(epoch * 1000);
    return date.getFullYear() + "/" + zeropad(date.getMonth() + 1) + "/" + zeropad(date.getDate()) + " " + zeropad(date.getHours()) + ":" + zeropad(date.getMinutes()) + ":" + zeropad(date.getSeconds());
}