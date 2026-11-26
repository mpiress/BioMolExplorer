"""Width-adaptive selection rows for settings dialogs.

The rows follow the width of the container instead of a fixed column grid:
each data column takes a share of the available space and only the actions
column keeps a compact fixed width, so a table with many columns never
overflows the settings dialog horizontally.  Long values are truncated by the
cell they are rendered in and keep the full text in their tooltip.

A ``DataTable`` model is kept beside the rendered rows: it carries the column
labels and the per-row metadata the dialogs and the tests read.
"""
import flet as ft

BORDER = '#E2E8F0'
HEADING = '#F1F5F9'


class SelectedItemsTable:
    """Render only existing rows; keep a DataTable model for row metadata."""

    def __init__(self, columns, weights=None, actions_width=108):
        self.columns = tuple(columns)
        # One weight per data column; the last column holds the row actions.
        self.weights = tuple(weights) if weights else (1,)*(len(self.columns)-1)
        self.actions_width = actions_width
        self.table = ft.DataTable(
            columns=[ft.DataColumn(ft.Text(label)) for label in columns], rows=[],
            column_spacing=16, heading_row_color=HEADING,
            data_row_min_height=56, data_row_max_height=56,
            divider_thickness=.5, show_checkbox_column=False)
        self.list = ft.Column(spacing=0,tight=True)
        # No horizontal scroll: the cells share the width the dialog offers.
        self.control = ft.Container(self.list,visible=False,border=ft.Border.all(1,BORDER),
                                    border_radius=10,clip_behavior=ft.ClipBehavior.HARD_EDGE)

    def line(self, controls, *, heading=False, color=None):
        cells = [ft.Container(control,expand=weight,alignment=ft.Alignment.CENTER_LEFT)
                 for control,weight in zip(controls,self.weights)]
        cells.append(ft.Container(controls[len(self.weights)],width=self.actions_width,
                                  alignment=ft.Alignment.CENTER_LEFT))
        return ft.Container(ft.Row(cells,spacing=8,vertical_alignment=ft.CrossAxisAlignment.CENTER),
            padding=ft.Padding.symmetric(horizontal=10,vertical=8),bgcolor=HEADING if heading else color,
            border=None if heading else ft.Border.only(top=ft.BorderSide(.5,BORDER)))

    def set_rows(self, rows):
        self.table.rows = rows
        self.control.visible = bool(rows)
        if not rows:
            self.list.controls = []
            return
        header = [ft.Text(label,size=12,weight=ft.FontWeight.W_600,max_lines=1,
                          overflow=ft.TextOverflow.ELLIPSIS) for label in self.columns]
        self.list.controls = [self.line(header,heading=True)]+[
            self.line([cell.content for cell in row.cells],color=row.color) for row in rows]
