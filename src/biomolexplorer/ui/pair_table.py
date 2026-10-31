"""Shared presentation of scientific pairs in redocking and Amber Score."""
import flet as ft


class PairTable:
    def __init__(self,columns,page,min_width=680,selectable=False,pin_actions=False):
        self.page=page;self.min_width=min_width;self.pin_actions=pin_actions
        self.action_width=112
        self.table=ft.DataTable(columns=[ft.DataColumn(ft.Text(label)) for label in columns],rows=[],
            column_spacing=36,heading_row_color='#F1F5F9',data_row_min_height=72,
            data_row_max_height=72,divider_thickness=.5,show_checkbox_column=selectable)
        self.actions_table=None
        if pin_actions:
            self.table.columns=self.table.columns[:-1]
            self.actions_table=ft.DataTable(columns=[ft.DataColumn(ft.Text(columns[-1]))],rows=[],
                width=self.action_width,column_spacing=0,horizontal_margin=4,
                heading_row_color='#F1F5F9',data_row_min_height=72,data_row_max_height=72,divider_thickness=.5)
            self.control=ft.Row([ft.Container(ft.Row([self.table],scroll=ft.ScrollMode.AUTO),expand=True),
                ft.Container(self.actions_table,width=self.action_width)],spacing=0,
                vertical_alignment=ft.CrossAxisAlignment.START,on_size_change=self.resize)
        else:self.control=ft.Row([self.table],scroll=ft.ScrollMode.ALWAYS,on_size_change=self.resize)

    def set_rows(self,rows):
        if not self.pin_actions:self.table.rows=rows;return
        self.actions_table.rows=[ft.DataRow(cells=[row.cells[-1]],color=row.color) for row in rows]
        for row in rows:row.cells=row.cells[:-1]
        self.table.rows=rows

    def resize(self,event):
        width=max(self.min_width,round(event.width)-(self.action_width if self.pin_actions else 0))
        if self.table.width!=width:
            self.table.width=width
            self.page.update()
