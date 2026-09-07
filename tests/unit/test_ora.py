import pandas as pd

from clustbuster.core.enrichment import AllClusterOraResult
from clustbuster.plotting.ora import ora_heatmap_figure, ora_heatmap_height


def test_ora_heatmap_clusters_pathways_and_source_clusters() -> None:
    values = pd.DataFrame(
        {
            "cluster_id": ["a", "a", "b", "b"],
            "cluster": ["0", "0", "1", "1"],
            "rank": [1, 2, 1, 2],
            "term": ["T activation", "Cytokine signaling", "T activation", "Phagocytosis"],
            "adjusted_p_value": [0.001, 0.02, 0.2, 0.003],
            "odds_ratio": [8.0, 4.0, 1.2, 7.0],
            "combined_score": [40.0, 12.0, 1.0, 35.0],
            "overlap_genes": ["CD3D", "IL7R", "CD3D", "LYZ"],
        }
    )
    result = AllClusterOraResult(
        library="GO_Biological_Process_2025",
        source="Enrichr",
        markers=pd.DataFrame(),
        values=values,
        failures=(),
    )
    figure = ora_heatmap_figure(result, top_n_pathways=3)

    assert figure.data[0].type == "heatmap"
    assert len(figure.data) == 1
    assert figure.layout.height == ora_heatmap_height(3)
    assert figure.layout.xaxis.title.text == "Cluster"
    assert figure.layout.yaxis.title.text == "Pathway"


    relabeled = ora_heatmap_figure(result, labels={"a": "T cell", "b": "T cell"})
    assert len(relabeled.data[0].x) == 2
    assert list(relabeled.layout.xaxis.ticktext) == ["T cell", "T cell"]
    assert result.values["cluster"].tolist() == ["0", "0", "1", "1"]
