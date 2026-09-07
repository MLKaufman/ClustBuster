"""Display labels without changing cluster identities or analysis groupings."""
from clustbuster.models import Workspace


def annotation_labels(workspace: Workspace, mode: str) -> dict[str, str]:
    return {
        record.cluster_id.serialized: (
            record.annotation.strip() or record.cluster_id.display
            if mode == "annotation" else record.cluster_id.display
        )
        for record in workspace.annotations.records()
    }
