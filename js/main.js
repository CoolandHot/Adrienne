var p_adjust = 0.05;

const random4chars = () => {
    return [...Array(4)].map(() => 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'[Math.floor(Math.random() * 62)]).join('');
}

function readCSVFile(file, reverse = true) {
    return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = event => {
            const csvData = event.target.result;
            const rows = csvData.split(/\r?\n/);
            const headers = rows.shift().replaceAll('"', '').split(',');
            headers[0] = 'gene'; // the first header is empty
            const records = rows.map(row => {
                const values = row
                    .replaceAll('"', '')
                    .split(',')
                const record = {};
                headers.forEach((header, index) => {
                    if (reverse & (header == "log2FoldChange" | header == "stat")) {
                        record[header] = -parseFloat(values[index]);
                    } else {
                        record[header] = values[index];
                    }
                });
                return record;
            });

            const filteredRecords = records.filter(record => record['padj'] < p_adjust);
            resolve(filteredRecords);
        };
        reader.onerror = error => reject(error);
        reader.readAsText(file);
    });
}



// ===========================================================
// ************************** echarts ***********************
// ===========================================================

var d3Container = d3.select("#container");
var vennChart = venn.VennDiagram()
    .width(500)
    .height(500);
var tooltip = d3.select("body")
    .append("div")
    .attr("class", "venntooltip");

const map_csv = () => {
    return new Promise((resolve, reject) => {
        var myPromises = getInputValues()
            .filter(v => typeof (v['file']) != 'undefined')
            .map(async v => {
                let records = await readCSVFile(v['file'], v['reverse']);
                // only reserve the up/down-regulation
                records = records.filter(record => v['up_regulate'] == 'up' ? record['log2FoldChange'] > 0 : record['log2FoldChange'] < 0);
                // get the remaining gene list
                let geneList = records.map(record => record['gene']);
                return { "name": `${v['name']}(${v['up_regulate']})`, "value": geneList, "records": records }
            });
        Promise.all(myPromises).then((results) => {
            resolve(results);
        });
    });

}

// *********************************
// ******* make intersection *******
// *********************************
function intersectList(...arrays) {
    return arrays.reduce((a, b) => a.filter(c => b.includes(c)));
}

/**
 * This function returns an array of all combinations of `m` elements from `arr`.
 * The order of the elements in each combination doesn’t matter.
 * @param {Array} arr - the array.
 * @param {number} m - The number of elements to combine as a group.
 * @returns {Array} An array of dictionary objects. {combination: Array(m), indices: Array(m)}
 *  `combination`, which is an array containing a combination of `m` elements from the input array;
 *  `indices`, which is an array containing the indices of the elements in that combination.
 */

function getCombinations(arr, m) {
    let result = [];
    let combination = Array(m).fill(0);
    let indices = Array(m).fill(0);
    function makeCombinations(arr, m, start, index) {
        if (index === m) {
            result.push({ combination: combination.slice(), indices: indices.slice() });
            return;
        }
        for (let i = start; i <= arr.length - 1 && arr.length - i >= m - index; ++i) {
            combination[index] = arr[i];
            indices[index] = i;
            makeCombinations(arr, m, i + 1, index + 1);
        }
    }
    makeCombinations(arr, m, 0, 0);
    return result;
}

/**
 * This function returns an array of all combinations of intersection.
 * @param {Array} dataset - the data table of records
 * @returns {Array} An array of dictionary objects. {sets: Array, size: number, records: Array}
 */
const compute_intersect = (dataset) => {
    // one set specially include the records
    let result = dataset.map(v => ({ sets: [v['name']], size: v['value'].length, records: v['records'] }))
    // two or more sets for intersections, the records are the intersect gene list, no detailed log2FC nor p.adjust
    for (let i = 1; i < dataset.length; i++) {
        let comb_res = getCombinations(dataset, i + 1)
            .map(v => {
                let intersect_v = intersectList(...v['combination'].map(x => x['value']))
                if (intersect_v.length > 0) {
                    return {
                        sets: v['combination'].map(x => x['name']),
                        size: intersect_v.length,
                        records: intersect_v
                    }

                }
            })
            .filter(x => typeof x != 'undefined')

        if (comb_res.length > 0) {
            result = result.concat(comb_res)
        }
    }
    return result
}


const plot_update = (dataset) => {
    document.getElementById("chart-title").innerText = dataset.map(v => v['name']).join(" - ");

    var data = compute_intersect(dataset);
    d3Container.datum(data).call(vennChart);

    // outline of each circle
    d3Container.selectAll("path")
        .style("stroke-opacity", 0)
        .style("stroke", "#fff")
        .style("stroke-width", 2);

    // add listeners to all the groups to display tooltip on mouseover
    d3Container.selectAll("g")
        .on("mouseover", function (d, i) {
            // sort all the areas relative to the current item
            venn.sortAreas(d3Container, d);
            // Display a tooltip with the current size
            tooltip.transition().duration(100).style("opacity", .9);
            tooltip.text(d.size + " genes");
            // highlight the current path
            var selection = d3.select(this).transition("tooltip").duration(200);
            selection.select("path")
                .style("stroke-width", 3)
                .style("fill-opacity", d.sets.length == 1 ? .4 : .1)
                .style("stroke-opacity", 1);
        })

        .on('click', (d, i) => {
            // sort all the areas relative to the current item
            venn.sortAreas(d3Container, d);
            // pop up detailed tables
            document.getElementById("panel").classList.remove("hidden");

            var records = d.records;
            var tableNames = d.sets, table;
            if (tableNames.length == 1) {
                table = '<h3>' + tableNames + '</h3>';
                table += records.length + ' genes<br>';
                table += '<table><thead><tr><th>gene</th><th>log2FoldChange</th><th>p.adjust</th></tr></thead><tbody>';
                for (let j = 0; j < records.length; j++) {
                    table += '<tr>' + `<td>${records[j]['gene']}</td>` + `<td>${records[j]['log2FoldChange']}</td>` + `<td>${records[j]['padj']}</td>` + `</tr>`;
                }
                table += '</tbody></table>'
            } else {
                table = '<h3>' + tableNames.join(" - ") + '</h3>';
                table += records.length + ' genes<br>';
                table += records.join(", ")
            }
            document.querySelector('.panel-content').innerHTML = table;

        })

        .on("mousemove", function () {
            tooltip.style("left", (d3.event.pageX) + "px")
                .style("top", (d3.event.pageY - 28) + "px");
        })

        .on("mouseout", function (d, i) {
            tooltip.transition().duration(200).style("opacity", 0);
            var selection = d3.select(this).transition("tooltip").duration(200);
            selection.select("path")
                .style("stroke-width", 0)
                .style("fill-opacity", d.sets.length == 1 ? .25 : .0)
                .style("stroke-opacity", 0);
        });

}

// ===========================================================
// **************** set information box control **************
// ===========================================================
function addBox() {
    let newNode = document
        .querySelector('.set-information')
        .cloneNode(true);
    // clear all information
    newNode.querySelector("input[type='file']").value = "";
    newNode.querySelector("input[type='text']").value = "";
    newNode.querySelector("span").innerText = "Choose File";
    newNode.querySelector(".dataFilter-paras input[type='checkbox']").checked = true;
    let new_radio_group = `regulation_${random4chars()}`;
    newNode.querySelectorAll(".dataFilter-paras input[type='radio']")
        .forEach(o => o.name = new_radio_group)

    document
        .querySelector('.web-gadgets')
        .appendChild(newNode);

    overwrite_listers();
}

function removeBox() {
    let parent_box = document.querySelector('.web-gadgets');
    if (parent_box.childElementCount > 3) {
        parent_box.removeChild(parent_box.lastChild);
    }
}


// ===========================================================
// *********************** listeners *******************
// ===========================================================
const overwrite_listers = () => {

    // file picker 
    document.querySelectorAll("input[type='file']")
        .forEach(dom => {
            dom.addEventListener('change', (event) => {
                const file = event.target.files[0];
                // fill in the name in the input
                let input_dom = event.target
                    .parentNode
                    .parentNode
                    .querySelector('input[type="text"][name="set-name"]');
                if (input_dom.value == "") {
                    input_dom.value = file.name.replace(".csv", "");
                    input_dom.dispatchEvent(new Event("change"));
                }
                // change the "Choose File" to the file name
                let span_dom = event.target
                    .parentNode
                    .querySelector('span');
                span_dom.innerText = file.name.slice(-20)
            });
        })

    // clear file picker and input name
    document.querySelectorAll('.clear-file')
        .forEach(i => i.addEventListener('click', event => {
            let parent = event.target.parentNode;
            parent.querySelectorAll("input")
                .forEach(input => input.value = "")
            parent.querySelector("span").innerText = "Choose File"
        }))

}

overwrite_listers();

const getInputValues = () => {
    let results = [];
    document.querySelectorAll(".set-information")
        .forEach(one => {
            let res = {};
            res['file'] = one.querySelector("input[type='file']").files[0];
            res['name'] = one.querySelector("input[type='text']").value;
            res['reverse'] = one.querySelector(".dataFilter-paras input[type='checkbox']").checked;
            res['up_regulate'] = one.querySelector(".dataFilter-paras input[type='radio']:checked").value;
            results.push(res)
        })
    return results
}

document.querySelector("#p_adjust")
    .addEventListener('change', event => {
        p_adjust = parseFloat(event.target.value)
    })


document.getElementById("panel").addEventListener('click', function (event) {
    if (!document.querySelector('.panel-content').contains(event.target)) {
        document.getElementById("panel").classList.add("hidden");
    }
});

document.getElementById("createChart").addEventListener('click', async event => {
    let data = await map_csv();
    plot_update(data)
});