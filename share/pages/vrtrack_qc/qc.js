var getQueryStringParameterByName = function(name) {
    name = name.replace(/[\[]/, "\\[").replace(/[\]]/, "\\]");
    var regex = new RegExp("[\\?&]" + name + "=([^&#]*)");
    var results = regex.exec(location.search);
    return results == null ? "" : decodeURIComponent(results[1].replace(/\+/g, " "));
};

// function to fill an array with the properties of a given label
var fillProperties = function(labelProperties, label, list) {
    var props = labelProperties[label];
    var arr = [];
    var i;
    var plen = props.length;
    for (i = 0; i < plen; i += 1) {
        arr.push(props[i]);
    }
    list(arr);
};

// lane information we get from getQCGraphData(nodes_of_label) needs to be
// massaged
var massageLaneProperties = function(lane, isGroupAdmin, aqs) {
    lane['new_qcgrind_qc_status'] = ko.observable(lane['qcgrind_qc_status']);
    lane['display_graphs'] = ko.observable(false);
    lane['display_stats'] = ko.observable(false);
    
    if (isGroupAdmin && ! lane.is_admin) {
        isGroupAdmin = false;
    }
    
    // create an auto_qc property that computes the
    // result so that the result can auto-update when
    // the user changes the settings. We need to put
    lane['auto_qc'] = ko.pureComputed(function() {
        var fails = [];
        
        // check we have data and end early if not
        if (! lane.lanelet_stats['alignmentstats:total length'] && ! lane.lanelet_stats['alignmentstats:sequences']) {
            fails.push(['Missing data', 'The alignment stats for this lanelet are missing.' ]);
        }
        else if (lane.lanelet_stats['alignmentstats:total length'] == 0 || lane.lanelet_stats['alignmentstats:sequences'] == 0) {
            fails.push(['Empty bam file check', 'The cram file provided for this lanelet contains no sequences.' ]);
        }
        else {
            // follow npg
            if (lane['NPG Status'] === 'fail') {
                fails.push(['NPG QC status check', 'The lane failed the NPG QC check, so we auto-fail as well since this data will not be auto-submitted to EGA/ENA.' ]);
            }
            
            // genotype check
            var gtype_regex_str = aqs['GType regex'].value();
            var gtype_regex = new RegExp(gtype_regex_str, "i");
            if (! gtype_regex.test(lane.genotype_status)) {
                fails.push(['Genotype check', 'The status (' + lane.genotype_status + ') does not match the regex /' + gtype_regex_str + '/.']);
            }
            var gtype_common_snp_count = lane['gtcheckdata:common_snp_count'];
            if (gtype_common_snp_count) {
                var gtype_min_s = parseInt(aqs['GType min sites'].value());
                if (parseInt(gtype_common_snp_count) < gtype_min_s) {
                    fails.push(['Genotype check', 'The number of common sites called (' + gtype_common_snp_count + ') is less than ' + gtype_min_s + '.']);
                }
            }
            var gtype_ratio_str = lane['gtcheckdata:ratio'];
            //*** should we have an auto qc setting for min ratio?
            if (gtype_common_snp_count || gtype_ratio_str) {
                var gtype_concordance_str = lane['gtcheckdata:concordance'];
                if (gtype_concordance_str) {
                    var gtype_min_c = parseFloat(aqs['GType min concordance'].value());
                    var gtype_concordance = parseFloat(gtype_concordance_str);
                    if (gtype_concordance < gtype_min_c) {
                        fails.push(['Genotype check', 'The concordance (' + gtype_concordance + ') is less than ' + gtype_min_c + '.']);
                    }
                }
            }
            
            // mapped bases
            var min_mapped_base_percentage = parseFloat(aqs['Mapped base %'].value());
            var clip_bases = lane.lanelet_stats['alignmentstats:total length'] - lane.lanelet_stats['alignmentstats:bases trimmed'];
            var bases_mapped_c = lane.lanelet_stats['alignmentstats:bases mapped (cigar)'];
            var mbtest = (100 * bases_mapped_c / clip_bases).toFixed(2);
            if (mbtest < min_mapped_base_percentage) {
                fails.push(['Mapped bases', 'Less than ' + min_mapped_base_percentage + '% bases mapped after clipping (' + mbtest + '%).']);
            }
            
            // duplicate reads
            var max_duplicate_read_percentage = parseFloat(aqs['Duplicate read %'].value());
            var reads_mapped = lane.lanelet_stats['alignmentstats:reads mapped'];
            if (reads_mapped) {
                // var reads_mapped_after_rmdup = lane['alignmentstats:reads mapped after rmdup'];
                // var dup_reads = reads_mapped - reads_mapped_after_rmdup;
                // var drtest = (100 * dup_reads / reads_mapped).toFixed(2);
                if (parseFloat(lane.lanelet_stats['Dup%']) > max_duplicate_read_percentage) {
                    fails.push(['Duplicate reads', 'More than ' + max_duplicate_read_percentage + '% reads were duplicates (' + lane.lanelet_stats['Dup%'] + '%).']);
                }
            }
            
            // properly paired mapped reads
            var min_mapped_reads_properly_paired_percentage = parseFloat(aqs['Mapped reads properly paired %'].value());
            if (reads_mapped) {
                var properly_paired = lane.lanelet_stats['alignmentstats:reads properly paired'];
                var pptest = (100 * properly_paired / reads_mapped).toFixed(2);
                if (pptest < min_mapped_reads_properly_paired_percentage) {
                    fails.push(['Reads mapped in a proper pair', 'Less than ' + min_mapped_reads_properly_paired_percentage + '% of reads that were mapped are in a proper pair (' + pptest + '%).']);
                }
            }
            
            // error rate
            var max_error_rate = parseFloat(aqs['Error rate'].value());
            var error_rate = parseFloat(lane.lanelet_stats['alignmentstats:error rate']); // or should this be the error rate % stored in lane['Error Rate']?
            if (error_rate && error_rate > max_error_rate) {
                fails.push(['Error rate', 'The error rate is higher than ' + max_error_rate + ' (' + +error_rate + ').']);
            }
            
            // number of insertions vs deletions
            var inum = lane.lanelet_stats['alignmentstats:number of insertions'];
            var dnum = lane.lanelet_stats['alignmentstats:number of deletions'];
            if (inum != null && dnum != null) {
                var max_ins_to_del_ratio = parseFloat(aqs['Max in/del ratio'].value());
                var min_ins_to_del_ratio = parseFloat(aqs['Min in/del ratio'].value());
                var ratio = dnum ? (inum / dnum).toFixed(2) : 0;
                if (!dnum || (ratio > max_ins_to_del_ratio)) {
                    fails.push(['InDel ratio', 'The Ins/Del ratio is bigger than ' + max_ins_to_del_ratio + ' (' + ratio + ').']);
                    $reason = "The Ins/Del ratio is bigger than $max ($inum/$dnum).";
                }
                else if (dnum && (!inum || (ratio < min_ins_to_del_ratio))) {
                    fails.push(['InDel ratio', 'The Ins/Del ratio is smaller than ' + min_ins_to_del_ratio + ' (' + ratio + ').']);
                }
            }
            
            // insert size
            if (lane['is_paired_read']) {
                if (parseInt(lane.lanelet_stats['alignmentstats:reads paired']) == 0) {
                    fails.push(['Insert size', 'Zero paired reads, yet flagged as paired']);
                }
                else if (parseFloat(lane.lanelet_stats['alignmentstats:insert size average']) == 0) {
                    fails.push(['Insert size', 'The insert size not available, yet flagged as paired']);
                }
                else {
                    // npg_cram_stats_parser step
                    // calculates these from all the
                    // (1000s of) IS values in the
                    // stats file, but only for fixed
                    // 80/25 settings. Sadly this causes
                    // a fail for everything, and since
                    // it's not user-configurable,
                    // excluding for now...
                    
                    // check how wide the insert size distribution is and
                    // whether 80% of the data lies within 25% from the max peak
                    // (e.g. [mpeak*(1-0.25),mpeak*(1+0.25)])
                    // only libraries can be failed based on wrong insert size. The
                    // lanes are always passed as long as the insert size is
                    // consistent with other lanes from the same library.
                    // var amount = lane.lanelet_stats['alignmentstats:inserts within 25% max peak'];
                    // var range  = lane.lanelet_stats['alignmentstats:peak window containing 80% of inserts'];
                    // if (amount != null && range != null) {
                    //     amount = parseFloat(amount).toFixed(2);
                    //     range = parseFloat(range).toFixed(2);
                    //     if (amount < 80) {
                    //         fails.push(['Insert size', 'Fail library, less than 80% of the inserts are within 25% of max peak (' + amount + '%).']);
                    //     }
                    //     if (range > 25) {
                    //         fails.push(['Insert size', 'Fail library, 80% of inserts are not within 25% of the max peak (' + range + '%).']);
                    //     }
                    // }
                    // else {
                    //     fails.push(['Insert size', 'The insert size not available, yet flagged as paired']);
                    // }
                }
            }
            
            // overlapping base duplicate percent
            // calculate the proportion of mapped bases duplicated e.g. if a fragment
            // is 160bp - then 40bp out of 200bp sequenced (or 20% of bases sequenced
            // in the fragment are duplicate sequence)
            //
            //------------->
            //          <------------
            //        160bp
            //|---------------------|
            //          |--|
            //          40bp
            var max_overlapping_base_duplicate_percent = parseFloat(aqs['Overlapping Base Duplicate %'].value());
            var dupl_mapped_bases = lane.lanelet_stats['alignmentstats:dupl mapped bases'];
            var total_mapped_bases = lane.lanelet_stats['alignmentstats:total mapped bases'];
            if (dupl_mapped_bases && total_mapped_bases) {
                var obtest = ((dupl_mapped_bases * 100) / total_mapped_bases).toFixed(2); // (same as lane['Overlap dup%'])
                if (obtest > max_overlapping_base_duplicate_percent) {
                    fails.push(['Overlap duplicate base percent', 'The percent of bases duplicated due to reads of a pair overlapping (' + obtest + '%) is greater than ' + max_overlapping_base_duplicate_percent + '%.']);
                }
            }
            
            // maximimum indels per cycle
            var max_ic_above_median = parseInt(aqs['Max indel/cycle above median'].value());
            var median_oicc = lane.lanelet_stats['alignmentstats:median of indel cycle counts'];
            var max_oicc = lane.lanelet_stats['alignmentstats:max of indel cycle counts'];
            if (median_oicc && max_oicc) {
                median_oicc = median_oicc.split(',');
                max_oicc = max_oicc.split(',');
                for (i = 0; i < 4; i += 1) {
                    if (max_oicc[i] > max_ic_above_median * median_oicc[i]) {
                        fails.push(['InDels per Cycle', 'Some ' + indelsPerCycleLookup[i] + ' per cycle exceed ' + max_ic_above_median + ' x median (max ' + max_oicc[i] + ' vs median ' + median_oicc[i] + ').']);
                    }
                }
            }
        }
        
        var status = 'pass';
        if (fails.length > 0) {
            status = 'fail';
        }
        return [status, fails];
    });
};

// function to call qc methods to get results from the graph db
var indelsPerCycleLookup = ['Insertion(fwd)', 'Insertion(rev)', 'Deletion(fwd)', 'Deletion(rev)'];
var qcStates = { 'Pass': 'passed', 'Fail': 'failed', 'Invst': 'investigate', 'GTPend': 'gt_pending', 'Pend': 'pending' };
var getQCGraphData = function(method, args, subargs, loading, errors) {
    loading.removeAll();
    errors.removeAll();
    
    vrpipeRestMethod('qc', method, args, loading, errors, function(data) {
        var resultStore = subargs.resultStore;
        // we won't push directly to this because even with
        // ratelimiting it's too slow; we'll push to temp arrays
        // and then set the resultStore to it
        var i;
        var len;
        var arr = [];
        var result;
        var type;
        
        switch (method) {
            case 'qc_website_admin':
                var adminUsers = subargs.adminUsers;
                var groupAdmins = subargs.groupAdmins;
                var gaArr = [];
                var groupConfig = subargs.groupConfig;
                var gcArr = [];
                len = data.length;
                for (i = 0; i < len; i += 1) {
                    result = data[i];
                    type = result.type;
                    delete result['type'];
                    switch (type) {
                        case 'admin':
                            if (adminUsers) {
                                adminUsers(result.users);
                            }
                            break;
                        case 'group_admins':
                            if (groupAdmins) {
                                gaArr.push(result);
                            }
                            break;
                        case 'group_config':
                            if (groupConfig) {
                                gcArr.push(result);
                            }
                            break;
                    }
                }
                
                if (groupAdmins) {
                    groupAdmins(gaArr);
                }
                if (groupConfig) {
                    groupConfig(gcArr);
                }
                break;
            
            case 'labels':
                resultStore.removeAll();
                var keys = [];
                Object.keys(data).forEach(function (key) {
                    keys.push(key);
                });
                keys.sort();
                var labelProperties = subargs.labelProperties;
                len = keys.length;
                var key;
                for (i = 0; i < len; i += 1) {
                    key = keys[i];
                    labelProperties[key] = data[key];
                    if (key !== 'Group' && key !== 'Study') {
                        arr.push(key);
                    }
                }      
                resultStore(arr);
                var viewLabels = subargs.viewLabels;
                viewLabels.push('Donor');
                viewLabels.push('Sample');
                viewLabels.push('Lanelet');
                
                // populate the groups select and the all properties
                fillProperties(labelProperties, 'Study', subargs.studyProperties);
                fillProperties(labelProperties, 'Donor', subargs.donorProperties);
                fillProperties(labelProperties, 'Sample', subargs.sampleProperties);
                getQCGraphData('nodes_of_label', { label: 'Group' }, { resultStore: subargs.groupNodes, sortProperty: 'name' }, loading, errors);
                break;
            
            case 'nodes_of_label':
                var flatten = subargs.flatten;
                var isGroupAdmin = false;
                if (args.label === 'Lanelet') {
                    isGroupAdmin = true;
                }
                var aqs = subargs.autoQCSettings;
                var slf = subargs.sampleLaneletsFilter;
                
                // (because we create computed functions for lanelets, we
                // can't use a normal for loop or they would all share the same
                // closure and things would be wrong)
                ko.utils.arrayForEach(data, function (datum) {
                    if (flatten) {
                        var props = datum.properties;
                        props.node_id = datum.id;
                        props.node_label = datum.label;
                        
                        if (datum.label === 'Lane') {
                            massageLaneProperties(props, isGroupAdmin, aqs);
                        }
                        else if (datum.label === 'Sample') {
                            props.display_lanelets = ko.observable(false);
                            
                            // massage the lanelets
                            var lanelets = [];
                            Object.keys(props.lanelet_nodes).forEach(function(key) {
                                var lanelet = this[key];
                                var laneletProps = lanelet.properties;
                                laneletProps.node_id = lanelet.id;
                                laneletProps.node_label = lanelet.label;
                                massageLaneProperties(laneletProps, isGroupAdmin, aqs);
                                lanelets.push(laneletProps);
                            }, props.lanelet_nodes);
                            props.lanelet_nodes = ko.observableArray(lanelets);
                            
                            // make the qc state numbers dynamic
                            Object.keys(qcStates).forEach(function(state) {
                                var status = qcStates[state];
                                props[state] = ko.pureComputed(function() {
                                    var n = 0;
                                    ko.utils.arrayForEach(props.lanelet_nodes(), function (lanelet) {
                                        if (lanelet.new_qcgrind_qc_status() == status) {
                                            n += 1;
                                        }
                                    });
                                    return n;
                                });
                            });
                            
                            // make dynamic properties that summarise the
                            // combined properties of the lanelets, depending
                            // on which lanelets pass a new_qcgrind_qc_status
                            // filter
                            props.lanelet_stats = ko.pureComputed(function() {
                                var stats = { 'Raw Gb': 0,  'Mapped Gb': 0, 'Mapped-dups Gb': 0, 'Net Gb': 0, 'alignmentstats:reads mapped': 0, 'alignmentstats:reads mapped after rmdup': 0, 'alignmentstats:bases mapped after rmdup': 0, 'alignmentstats:dupl mapped bases': 0, 'alignmentstats:total mapped bases': 0, 'alignmentstats:bases mapped (cigar)': 0, 'alignmentstats:mismatches': 0 };
                                var i;
                                var laneletsLength = lanelets.length;
                                var lanelet;
                                var ls;
                                for (i = 0; i < laneletsLength; i += 1) {
                                    lanelet = lanelets[i];
                                    if (slf[lanelet.new_qcgrind_qc_status()]()) {
                                        continue;
                                    }
                                    ls = lanelet.lanelet_stats;
                                    
                                    stats['Raw Gb']         += ls['Raw Gb'] != -1 ? +ls['Raw Gb'] : 0;
                                    stats['Mapped Gb']      += ls['Mapped Gb'] != -1 ? +ls['Mapped Gb'] : 0;
                                    stats['Mapped-dups Gb'] += ls['Mapped-dups Gb'] != -1 ? +ls['Mapped-dups Gb'] : 0;
                                    stats['Net Gb']         += ls['Net Gb'] != -1 ? +ls['Net Gb'] : 0;
                                    
                                    // to calculate Dup%, Overlap dup%, Depth
                                    // and error rate later, we need the
                                    // following:
                                    stats['alignmentstats:reads mapped']             += +ls['alignmentstats:reads mapped'];
                                    stats['alignmentstats:reads mapped after rmdup'] += +ls['alignmentstats:reads mapped after rmdup'];
                                    stats['alignmentstats:bases mapped after rmdup'] += +ls['alignmentstats:bases mapped after rmdup'];
                                    stats['alignmentstats:dupl mapped bases']        += +ls['alignmentstats:dupl mapped bases'];
                                    stats['alignmentstats:total mapped bases']       += +ls['alignmentstats:total mapped bases'];
                                    stats['alignmentstats:bases mapped (cigar)']     += +ls['alignmentstats:bases mapped (cigar)'];
                                    stats['alignmentstats:mismatches']               += +ls['alignmentstats:mismatches'];
                                }
                                
                                // finalise the summary stats now that we have
                                // all the filtered stats
                                stats['Raw Gb'] = stats['Raw Gb'].toFixed(2);
                                stats['Mapped Gb'] = stats['Mapped Gb'].toFixed(2);
                                stats['Mapped-dups Gb'] = stats['Mapped-dups Gb'].toFixed(2);
                                stats['Net Gb'] = stats['Net Gb'].toFixed(2);
                                
                                var reads_mapped             = stats['alignmentstats:reads mapped'];
                                var reads_mapped_after_rmdp  = stats['alignmentstats:reads mapped after rmdup'];
                                var bases_mapped_after_rmdup = stats['alignmentstats:bases mapped after rmdup'];
                                var dupl_mapped_bases        = stats['alignmentstats:dupl mapped bases'];
                                var total_mapped_bases       = stats['alignmentstats:total mapped bases'];
                                var bases_mapped_cigar       = stats['alignmentstats:bases mapped (cigar)'];
                                if (reads_mapped && reads_mapped_after_rmdp) {
                                    stats['Dup%'] = rounder((100 / reads_mapped) * (reads_mapped - reads_mapped_after_rmdp));
                                }
                                else {
                                    stats['Dup%'] = -1;
                                }
                                if (dupl_mapped_bases && total_mapped_bases) {
                                    stats['Overlap dup%'] = rounder((dupl_mapped_bases * 100) / total_mapped_bases);
                                }
                                else {
                                    stats['Overlap dup%'] = -1;
                                }
                                if (bases_mapped_cigar) {
                                    stats['Error Rate'] = rounder((100 / bases_mapped_cigar) * stats['alignmentstats:mismatches']);
                                }
                                else {
                                    stats['Error Rate'] = -1;
                                }
                                if (bases_mapped_after_rmdup) {
                                    // *** qcgrind sample_view.pl does something odd
                                    // where the mapping raw_bases is expected to be
                                    // different from the lane raw_bases; is this an
                                    // exome thing??
                                    stats['Depth'] = rounder(bases_mapped_after_rmdup / 3000000000); // human hard-coding
                                }
                                else {
                                    stats['Depth'] = -1;
                                }
                                
                                return stats;
                            });
                        }
                        
                        arr.push(props);
                    }
                    else {
                        arr.push(datum);
                    }
                });
                resultStore(arr);
                if (subargs.isGroupAdmin) {
                    subargs.isGroupAdmin(isGroupAdmin);
                }
                break;
            
            case 'node_by_id':
                result = subargs.result;
                result(data);
                break;
            
            case 'donor_qc':
                var donorAdmin = subargs.donorAdmin;
                var donorSampleStatus = subargs.donorSampleStatus;
                var donorGenderResults = subargs.donorGenderResults;
                var donorFluidigmDiscordance = subargs.donorFluidigmDiscordance;
                var donorGenotypingDiscordance = subargs.donorGenotypingDiscordance;
                var donorCopyNumberSummary = subargs.donorCopyNumberSummary;
                var donorAberrantRegions = subargs.donorAberrantRegions;
                var donorAberrantPolysomy = subargs.donorAberrantPolysomy;
                var donorCopyNumberPlot = subargs.donorCopyNumberPlot;
                var donorLOHCalls = subargs.donorLOHCalls;
                var donorPluritestSummary = subargs.donorPluritestSummary;
                var donorPluritestPlots = subargs.donorPluritestPlots;
                donorSampleStatus.removeAll();
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
                var ssArr = [];
                var genderArr = [];
                var discFArr = [];
                var discGArr = [];
                var cnsArr = [];
                var arArr = [];
                var apArr = [];
                var lohArr = [];
                var psArr = [];
                var ppArr = [];
                
                len = data.length;
                for (i = 0; i < len; i += 1) {
                    result = data[i];
                    type = result.type;
                    delete result['type'];
                    switch (type) {
                        case 'admin':
                            donorAdmin(result);
                            break;
                        case 'sample_status':
                            result['new_qc_failed_reason'] = ko.observable(result['qc_failed_reason']);
                            result['new_qc_status'] = ko.observable(result['qc_status']);
                            result['new_qc_passed_fluidigm'] = ko.observable(result['qc_passed_fluidigm']);
                            result['new_qc_exclude_from_analysis'] = ko.observable(result['qc_exclude_from_analysis']);
                            ssArr.push(result);
                            break;
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
                            if (result.hasOwnProperty('graph')) {
                                result['graph'] = '/file' + result.graph;
                                if (result.graph === '/filenull') {
                                    result['graph'] = '-none-';
                                }
                            }
                            arArr.push(result);
                            break;
                        case 'aberrant_polysomy':
                            if (result.hasOwnProperty('graph')) {
                                result['graph'] = '/file' + result.graph;
                            }
                            apArr.push(result);
                            break;
                        case 'copy_number_plot':
                            donorCopyNumberPlot('/file' + result.plot);
                            break;
                        case 'loh_calls' :
                            lohArr.push(result);
                            break;
                        case 'pluritest_summary' :
                            psArr.push(result);
                            break;
                        case 'pluritest_plot':
                            result['path'] = '/file' + result.path;
                            ppArr.push(result);
                            break;
                    }
                }
                
                donorSampleStatus(ssArr);
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
                var sdResults = subargs.sdResults;
                sdResults.removeAll();
                len = data.length;
                for (i = 0; i < len; i += 1) {
                    arr.push(data[i]);
                }
                sdResults(arr);
                break;
            
            default:
                errors.push('invalid qc method: ' + method);
        }
    }, subargs);
};