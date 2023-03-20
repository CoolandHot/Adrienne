var file_paras = [
    { "file": "", "name": "", "reverse": true },
    { "file": "", "name": "", "reverse": true },
    { "file": "", "name": "", "reverse": true },
    { "file": "", "name": "", "reverse": true },
];
var p_adjust = 0.05;


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
// ***************** file picker listeners ******************
// ===========================================================
const filePickers = ['filePicker1', 'filePicker2', 'filePicker3', 'filePicker4'];
for (const picker of filePickers) {
    document.getElementById(picker).addEventListener('click', () => {
        const input = document.createElement('input');
        input.type = 'file';
        input.accept = '.csv';
        input.onchange = event => {
            const file = event.target.files[0];
            let num = parseInt(picker.replace("filePicker", ""));
            let input_dom = document.getElementById(`file${num}`);
            if (input_dom.value == "") {
                input_dom.value = file.name.replace(".csv", "");
            }
            file_paras[num - 1] = { "file": file, "name": input_dom.value, "reverse": document.getElementById(`reverse_fc${num}`).checked }
        };
        input.click();
    });
}
// ===========================================================
// *********************** input listeners *******************
// ===========================================================
const inputs = document.querySelectorAll('input');
for (const input of inputs) {
    input.addEventListener('change', event => {
        if (event.target.type == "text") {
            if (event.target.id == "p_adjust") {
                p_adjust = parseFloat(document.getElementById("p_adjust").value)
            }
            if (event.target.id.includes('file')) { // file name
                file_paras[parseInt(event.target.id.replace("file", "")) - 1]["name"] = event.target.value
            }
        }
        if (event.target.type == "checkbox")
            file_paras[parseInt(event.target.id.replace("reverse_fc", "")) - 1]["reverse"] = event.target.checked;
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


// var chartOption = {
//     tooltip: {
//         trigger: 'item',
//         formatter: function (params) {
//             var records = params.data.records;
//             var table = '<h3>' + params.marker + ' ' + params.name + '</h3>';
//             table += '<table><thead><tr><th>gene</th><th>log2FoldChange</th><th>p.adjust</th></tr></thead><tbody>';
//             for (let j = 0; j < records.length; j++) {
//                 table += '<tr>' + `<td>${records[j]['gene']}</td>` + `<td>${records[j]['log2FoldChange']}</td>` + `<td>${records[j]['padj']}</td>` + `</tr>`;
//             }
//             table += '</tbody></table>';
//             return table;
//         }
//     },
//     toolbox: {
//         show: true,
//         feature: {
//             mark: { show: true },
//             dataView: { show: true, readOnly: false },
//             restore: { show: true },
//             saveAsImage: { show: true }
//         }
//     },
//     calculable: false,
// };

const map_csv = (up_down = true) => {
    return new Promise((resolve, reject) => {
        var myPromises = file_paras
            .filter(v => v['file'] != "")
            .map(async v => {
                let records = await readCSVFile(v['file'], v['reverse']);
                // only reserve the up/down-regulation
                records = records.filter(record => up_down ? record['log2FoldChange'] > 0 : record['log2FoldChange'] < 0);
                // get the remaining gene list
                let geneList = records.map(record => record['gene']);
                return { "name": v['name'], "value": geneList, "records": records }
            });
        Promise.all(myPromises).then((results) => {
            resolve(results);
        });
    });

}

const plot_update = (dataset, up_down) => {
    document.getElementById("chart-title").innerText = dataset.map(v => v['name']).join(" - ") + up_down ? "up-regulated" : "down-regulated";

    // *********************************
    // ******* make intersection *******
    // *********************************
    // singleton
    let data_one = dataset.map(v => ({ sets: [v['name']], size: v['value'].length, records: v['records'] }));

    // intersection of two sets
    let data_two01 = dataset[0]['value'].filter(item => dataset[1]['value'].includes(item)), data_two02, data_two12, data_two03, data_two13, data_two23, data_two = [
        {
            sets: [dataset[0]['name'], dataset[1]['name']],
            size: data_two01.length,
            records: data_two01
        }
    ];

    // intersection of three/four sets
    let data_three = [], data_four = [];

    if (dataset.length > 2) {
        data_two02 = dataset[0]['value'].filter(item => dataset[2]['value'].includes(item));
        data_two12 = dataset[1]['value'].filter(item => dataset[2]['value'].includes(item));
        data_two.concat([
            { sets: [dataset[0]['name'], dataset[2]['name']], size: data_two02.length, records: data_two02 },
            { sets: [dataset[1]['name'], dataset[2]['name']], size: data_two12.length, records: data_two12 }
        ])
        let data_three012 = dataset[0]['value'].filter(item => dataset[1]['value'].includes(item) && dataset[2]['value'].includes(item));
        data_three.concat([
            {
                sets: dataset.slice(0, 3).map(v => v['name']),
                size: data_three012.length,
                records: data_three012
            }
        ])
    }
    if (dataset.length > 3) {
        data_two03 = dataset[0]['value'].filter(item => dataset[3]['value'].includes(item));
        data_two13 = dataset[1]['value'].filter(item => dataset[3]['value'].includes(item));
        data_two23 = dataset[2]['value'].filter(item => dataset[3]['value'].includes(item));
        data_two.concat([
            { sets: [dataset[0]['name'], dataset[3]['name']], size: data_two03.length, records: data_two03 },
            { sets: [dataset[1]['name'], dataset[3]['name']], size: data_two13.length, records: data_two13 },
            { sets: [dataset[2]['name'], dataset[3]['name']], size: data_two23.length, records: data_two23 }
        ])
        let data_three013 = dataset[0]['value'].filter(item => dataset[1]['value'].includes(item) && dataset[3]['value'].includes(item));
        let data_three023 = dataset[0]['value'].filter(item => dataset[2]['value'].includes(item) && dataset[3]['value'].includes(item));
        let data_three123 = dataset[1]['value'].filter(item => dataset[2]['value'].includes(item) && dataset[3]['value'].includes(item));
        data_three.concat([
            {
                sets: [0, 1, 3].map(ind => dataset[ind]).map(v => v['name']),
                size: data_three013.length,
                records: data_three013
            },
            {
                sets: [0, 2, 3].map(ind => dataset[ind]).map(v => v['name']),
                size: data_three023.length,
                records: data_three023
            },
            {
                sets: [1, 2, 3].map(ind => dataset[ind]).map(v => v['name']),
                size: data_three123.length,
                records: data_three123
            },
        ])

        let data_four0123 = dataset[0]['value'].filter(item => dataset[1]['value'].includes(item) && dataset[2]['value'].includes(item) && dataset[3]['value'].includes(item));

        data_four.concat([
            {
                sets: dataset.map(v => v['name']),
                size: data_four0123.length,
                records: data_four0123
            }
        ])
    }
    // *********************************
    // *** end of make intersection ****
    // *********************************

    var data = data_one.concat(data_two).concat(data_three).concat(data_four);

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
            tooltip.text(d.size + " users");
            // highlight the current path
            var selection = d3.select(this).transition("tooltip").duration(200);
            selection.select("path")
                .style("stroke-width", 3)
                .style("fill-opacity", d.sets.length == 1 ? .4 : .1)
                .style("stroke-opacity", 1);
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



document.getElementById("down-regVisual").addEventListener('click', async event => {
    let data = await map_csv(false);
    plot_update(data, false)
});
document.getElementById("up-regVisual").addEventListener('click', async event => {
    let data = await map_csv(true);
    plot_update(data, true)
});