// Virüs konfigürasyonları
const configurations = {
    sars2: {
        name: "SARS-CoV-2",
        pdb: "data/7K8S.pdb",
        top3Json: "data/sars2_final_top3.json",
        structuralDesc: "<strong>7K8S (Spike Proteini)</strong> 3D analizi. Renkli bölgeler, tespit edilen farklı hedeflerdir.",
        typeMapping: { 0: "Aday Bölge 1", 1: "Aday Bölge 2", 2: "Aday Bölge 3" }
    },
    sars1: {
        name: "SARS-CoV-1",
        pdb: "data/5X58.pdb",
        top3Json: "data/sars1_final_top3.json",
        structuralDesc: "<strong>5X58</strong> 3D analizi. SARS-CoV-1 için evrimsel korunmuş bölgeler.",
        typeMapping: { 0: "Aday Bölge 1", 1: "Aday Bölge 2", 2: "Aday Bölge 3" }
    },
    flu: {
        name: "Influenza A",
        pdb: "data/1RVX.pdb",
        top3Json: "data/flu_final_top3.json",
        structuralDesc: "<strong>1RVX (Hemagglutinin)</strong> 3D analizi. Antijenik sürüklenmeye dirençli adaylar.",
        typeMapping: { 0: "Aday Bölge 1", 1: "Aday Bölge 2", 2: "Aday Bölge 3" }
    }
};

const candidateColorsRGB = [
    { r: 239, g: 68, b: 68 },  // Kırmızı
    { r: 59, g: 130, b: 246 },  // Mavi
    { r: 16, g: 185, b: 129 }   // Yeşil
];

let viewerInstance = null;
let currentTop3Data = [];
let currentVirusKey = 'sars2';

$(function () {
    viewerInstance = new PDBeMolstarPlugin();
    updateInterface('sars2');
});

async function updateInterface(virusKey) {
    currentVirusKey = virusKey;
    const config = configurations[virusKey];
    if (window.MathJax) MathJax.typesetPromise();

    const grid = document.getElementById('top-candidates-grid');
    grid.innerHTML = '<p>Yükleniyor...</p>';

    try {
        // Top3 JSON verisini oku
        const response = await fetch(config.top3Json);
        currentTop3Data = await response.json();

        // Ağırlıkları güncelle
        const w_s = currentTop3Data[0]?.w_stats ?? 0.5;
        const w_a = currentTop3Data[0]?.w_ai ?? 0.5;
        const formulaStr = "S_final = " + w_s.toFixed(2) + " * Stats + " + w_a.toFixed(2) + " * AI";
        document.getElementById('formula-display').innerText = formulaStr;
        document.getElementById('ml-status').innerText =
            "* " + config.name + " için otomatik hesaplanmış ağırlıklar: Stats=" + w_s.toFixed(2) + ", AI=" + w_a.toFixed(2);

        // Ground truth verilerini çek (IEDB epitopları)
        let groundTruthEpitopes = [];
        let fullGtData = {};
        try {
            const gtResponse = await fetch('data/ground_truth.json');
            fullGtData = await gtResponse.json();
            groundTruthEpitopes = fullGtData[virusKey] || [];
        } catch (e) {
            console.warn('Ground truth yüklenemedi:', e);
        }

        grid.innerHTML = '';

        currentTop3Data.forEach(function (candidate, index) {
            const rgb = candidateColorsRGB[index];
            const hexColor = "rgb(" + rgb.r + "," + rgb.g + "," + rgb.b + ")";
            const card = document.createElement('div');
            card.className = 'metric-card';
            card.style.borderTopColor = hexColor;

            const seqStr = candidate.seq_snippet || "N/A";
            const startPos = candidate.pos - 2;
            const endPos = candidate.pos + 2;

            card.innerHTML =
                '<span class="badge" style="background-color:' + hexColor + ';color:white;border:none;">Aday ' + (index + 1) + '</span>'
                + '<h4>' + config.typeMapping[index] + '</h4>'
                + '<div style="background:#f8fafc;padding:0.5rem;border-radius:4px;border:1px solid #e2e8f0;font-family:monospace;font-size:1.1rem;text-align:center;margin:10px 0;color:' + hexColor + ';font-weight:bold;letter-spacing:2px;">' + seqStr + '</div>'
                + '<p style="font-size:0.9rem;margin-bottom:0.5rem;">Pozisyon: <strong>' + startPos + ' – ' + endPos + '</strong> (Merkez: ' + candidate.pos + ')</p>'
                + '<p style="font-size:0.9rem;margin-bottom:0.25rem;">Erişilebilirlik (SASA): ' + candidate.sasa + ' Å²</p>'
                + '<p style="font-size:0.9rem;margin-bottom:0.25rem;">'
                + '<span title="ESM-2 protein dil modelinin bu pozisyonu ne kadar değiştirilemez bulduğu." style="cursor:help;border-bottom:1px dashed #94a3b8;">ESM-2 Kararlılık:</span>'
                + ' %' + (candidate.ai_score * 100).toFixed(1)
                + '</p>'
                + '<div class="metric-score" style="margin-top:1rem;">' + candidate.final_score.toFixed(3) + '<span>S<sub>final</sub></span></div>';

            grid.appendChild(card);
        });

        // ── Karşılaştırma panelini render et ──
        renderComparisonPanel(currentTop3Data, fullGtData, virusKey);

    } catch (err) {
        grid.innerHTML = '<p>Hata: ' + err.message + '</p>';
    }

    loadProteinStructure(config.pdb);
}

/**
 * Gelişmiş Sekans Eşleşme Mekanizması
 *
 * Sırasıyla şunları dener:
 * 1. Hedef virüste TAM eşleşme.
 * 2. Diğer virüslerde (örn. SARS-CoV-2) TAM eşleşme (Homoloji/Çapraz Reaktivite).
 * 3. Hedef virüste 1 amino asit farkla (Fuzzy) eşleşme.
 */
function calcBestSeqMatch(candidate, currentVirusKey, fullGtData) {
    const rawSeq = (candidate.seq_snippet || '').replace(/-/g, '').toUpperCase();
    if (!rawSeq || rawSeq.length === 0) return { score: 0, name: null, matchedEpitope: null, iedbId: null, matchType: 'none' };

    let bestMatch = { score: 0, name: null, matchedEpitope: null, iedbId: null, matchType: 'none', isConserved: false };

    // Tüm virüsleri tara (önce mevcut virüsü tara)
    const virusKeys = Object.keys(fullGtData);
    // Sıralama: Mevcut virüs en başta olsun
    const sortedKeys = [currentVirusKey, ...virusKeys.filter(k => k !== currentVirusKey)];

    for (const vKey of sortedKeys) {
        const epitopes = fullGtData[vKey] || [];

        for (const epitope of epitopes) {
            const epitopeSeq = (epitope.name || '').toUpperCase();
            if (!epitopeSeq || epitopeSeq.length < rawSeq.length) continue;

            let score = 0;
            let currentMatchType = 'none';

            // 1. Tam Substring Eşleşmesi
            if (epitopeSeq.indexOf(rawSeq) !== -1) {
                score = rawSeq.length / epitopeSeq.length;
                currentMatchType = (vKey === currentVirusKey) ? 'exact' : 'homology';
            }
            // 2. Fuzzy Eşleşme
            else if (score === 0) {
                const winLen = rawSeq.length;
                for (let i = 0; i <= epitopeSeq.length - winLen; i++) {
                    const window = epitopeSeq.slice(i, i + winLen);
                    let mismatches = 0;
                    for (let j = 0; j < winLen; j++) {
                        if (rawSeq[j] !== window[j]) mismatches++;
                    }
                    if (mismatches <= 1) { 
                        const fuzzyScore = (winLen - mismatches) / epitopeSeq.length;
                        if (fuzzyScore > score) {
                            score = fuzzyScore;
                            currentMatchType = (vKey === currentVirusKey) ? 'fuzzy' : 'homology-fuzzy';
                        }
                    }
                }
            }

            if (score > bestMatch.score) {
                bestMatch = {
                    score: score,
                    name: epitope.name,
                    matchedEpitope: epitopeSeq,
                    iedbId: epitope.iedb_id,
                    matchType: currentMatchType,
                    originVirus: vKey,
                    isConserved: false
                };
                if (currentMatchType === 'exact' && score > 0.8) break;
            }
        }
        if (bestMatch.matchType === 'exact' && bestMatch.score > 0.8) break;
    }

    // Konservasyon Kontrolü (SARS-CoV-1 vs SARS-CoV-2)
    if (bestMatch.score > 0 && (currentVirusKey === 'sars1' || currentVirusKey === 'sars2')) {
        const otherKey = (currentVirusKey === 'sars1') ? 'sars2' : 'sars1';
        const otherEpitopes = fullGtData[otherKey] || [];
        const isFoundInOther = otherEpitopes.some(e => (e.name || '').toUpperCase().indexOf(rawSeq) !== -1);
        if (isFoundInOther) {
            bestMatch.isConserved = true;
        }
    }

    return bestMatch;
}

function renderComparisonPanel(top3Data, fullGtData, currentVirusKey) {
    const compGrid = document.getElementById('comparison-grid');
    if (!compGrid) return;

    compGrid.innerHTML = '';

    top3Data.forEach(function (candidate, index) {
        const rgb = candidateColorsRGB[index];
        const hexColor = "rgb(" + rgb.r + "," + rgb.g + "," + rgb.b + ")";
        const result = calcBestSeqMatch(candidate, currentVirusKey, fullGtData);
        const pct = (result.score * 100).toFixed(1);
        const seqStr = candidate.seq_snippet || 'N/A';
        const rawSeq = seqStr.replace(/-/g, '');

        // Match Type Labels
        let label = '✗ IEDB\'de Bulunamadı';
        let barColor = '#94a3b8';
        let matchHint = '';

        if (result.matchType === 'exact' || result.matchType === 'homology') {
            label = '✓ Tam Epitop Eşleşmesi';
            barColor = result.score >= 0.5 ? '#10b981' : '#f59e0b';
            
            // SARS ailesi (SARS1/2) içindeki tüm tam eşleşmeleri 'korunmuş bölge' olarak işaretle
            if ((currentVirusKey === 'sars1' || currentVirusKey === 'sars2') && 
                (result.originVirus === 'sars1' || result.originVirus === 'sars2')) {
                matchHint = '*Bu bölge her iki virüs türünde de (SARS-CoV-1/2) bulunan, evrimsel olarak korunmuş bir epitoptur.';
            }
        } else if (result.matchType === 'fuzzy' || result.matchType === 'homology-fuzzy') {
            label = '⚠ Kısmi Eşleşme';
            barColor = '#f59e0b';
            matchHint = '*1 amino asit farkı ile eşleşme sağlandı.';
            if (result.matchType === 'homology-fuzzy') {
                matchHint += ' (' + result.originVirus.toUpperCase() + ' kaynağı)';
            }
        }

        const epitopeDisplay = result.name
            ? result.name.slice(0, 30) + (result.name.length > 30 ? '…' : '')
            : '—';

        // Eşleşen epitop içinde adayı highlight et
        let highlightedEpitope = '';
        if (result.matchedEpitope) {
            const idx = result.matchedEpitope.indexOf(rawSeq);
            if (idx !== -1) {
                highlightedEpitope =
                    result.matchedEpitope.slice(0, idx)
                    + '<mark style="background:#fef08a;padding:0 2px;border-radius:2px;font-weight:700;">'
                    + result.matchedEpitope.slice(idx, idx + rawSeq.length)
                    + '</mark>'
                    + result.matchedEpitope.slice(idx + rawSeq.length);
            } else {
                highlightedEpitope = result.matchedEpitope;
            }
        }

        const row = document.createElement('div');
        row.className = 'comparison-row';

        const iedbLink = result.iedbId
            ? '<a href="https://www.iedb.org/epitope/' + result.iedbId + '" target="_blank" style="text-decoration:none; color:#3b82f6; font-weight:700;">[IEDB:' + result.iedbId + ']</a>'
            : '';

        row.innerHTML =
            '<div class="comparison-candidate-label">'
            + '<span class="badge" style="background-color:' + hexColor + ';">Aday ' + (index + 1) + '</span>'
            + '<span class="seq-mono" style="color:' + hexColor + ';">' + seqStr + '</span>'
            + '<span style="font-size:0.7rem;color:#94a3b8;">(' + rawSeq + ')</span>'
            + '</div>'
            + '<div class="comparison-bar-wrap">'
            + '<div class="bar-label">'
            + '<span style="color:' + barColor + ';font-weight:600;">' + label + '</span>'
            + '<span class="epitope-name" style="font-family:monospace;">' + epitopeDisplay + ' ' + iedbLink + '</span>'
            + '</div>'
            + '<div class="comparison-bar-track">'
            + '<div class="comparison-bar-fill" style="width:0%;background:' + barColor + ';" data-target="' + pct + '"></div>'
            + '</div>'
            + (highlightedEpitope
                ? '<div style="margin-top:0.4rem;font-size:0.75rem;color:#64748b;font-family:monospace;">'
                + 'Epitop: ' + highlightedEpitope + ' ' + (matchHint ? '<br><span style="color:#3b82f6;font-size:0.7rem;">' + matchHint + '</span>' : '') + '</div>'
                : '<div style="margin-top:0.4rem;font-size:0.75rem;color:#94a3b8;font-style:italic;">Hiçbir IEDB epitobunda tam eşleşme yok.</div>')
            + '</div>'
            + '<div class="comparison-iou-badge" style="color:' + barColor + ';">'
            + pct + '%'
            + '<small>Epitop Kapsama</small>'
            + '</div>';

        compGrid.appendChild(row);
    });

    // Animasyon: bar'ları 0'dan hedef genişliğe aç
    setTimeout(function () {
        compGrid.querySelectorAll('.comparison-bar-fill').forEach(function (bar) {
            bar.style.width = bar.dataset.target + '%';
        });
    }, 100);
}


function loadProteinStructure(pdbPath) {
    const options = {
        customData: { url: pdbPath, format: 'pdb' },
        assemblyId: '1',
        bgColor: { r: 255, g: 255, b: 255 },
        hideControls: false,
        hideSequencePanel: false,
        sequencePanel: true,
        landscape: true,
        visualStyle: 'cartoon'
    };

    const viewerContainer = document.getElementById('myViewer');
    viewerInstance.render(viewerContainer, options);

    setTimeout(function () {
        if (viewerInstance && viewerInstance.plugin && viewerInstance.plugin.layout) {
            viewerInstance.plugin.layout.setProps({
                showSequence: true,
                regionState: { top: 'show', right: 'show' }
            });
        }
        renderMolstarStyles();
    }, 2000);

    if (viewerInstance.events && viewerInstance.events.loadComplete) {
        viewerInstance.events.loadComplete.subscribe(function () {
            setTimeout(renderMolstarStyles, 1000);
        });
    }
}

function renderMolstarStyles() {
    if (!viewerInstance || !viewerInstance.visual) return;

    const selectionData = [];
    currentTop3Data.forEach(function (candidate, index) {
        for (let i = -2; i <= 2; i++) {
            selectionData.push({
                struct_asym_id: 'A',
                residue_number: candidate.pos + i,
                color: candidateColorsRGB[index],
                focus: false
            });
        }
    });

    viewerInstance.visual.select({
        data: selectionData,
        nonSelectedColor: { r: 220, g: 220, b: 220 }
    });
}

function startAutomatedAnalysis() {
    const selected = document.getElementById('virus-select').value;
    const overlay = document.getElementById('loading-overlay');
    const fill = document.getElementById('progress-fill');
    const statusText = document.getElementById('status-text');

    overlay.style.display = 'flex';
    fill.style.width = '0%';

    const steps = [
        "Veri Setleri Yükleniyor...",
        "MSA Hizalaması Yapılıyor...",
        "ESM-2 Analizi Yürütülüyor...",
        "SASA ve Yüzey Erişilebilirliği Hesaplanıyor...",
        "IEDB Epitop Karşılaştırması Yapılıyor..."
    ];

    let step = 0;
    const interval = setInterval(function () {
        if (step < steps.length) {
            statusText.innerText = steps[step];
            fill.style.width = (((step + 1) / steps.length) * 100) + "%";
            step++;
        } else {
            clearInterval(interval);
            setTimeout(function () {
                overlay.style.display = 'none';
                updateInterface(selected);
            }, 500);
        }
    }, 800);
}

function downloadCurrentAnalysis() {
    const config = configurations[currentVirusKey];
    const filename = "smartepitope_" + currentVirusKey + "_top3.json";

    fetch(config.top3Json)
        .then(function (res) { return res.json(); })
        .then(function (data) {
            const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = filename;
            a.click();
            URL.revokeObjectURL(url);
        })
        .catch(function (err) { alert('İndirme hatası: ' + err.message); });
}
