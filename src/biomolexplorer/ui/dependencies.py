"""Keep dependent settings visible and preserve their values while inactive."""
import flet as ft


def bind_dependencies(controls, rules, writable, page):
    """Rules map child keys to predicates; chain existing parent callbacks."""
    def sync():
        for key, active in rules.items():
            if key in controls:
                controls[key].disabled = not writable or not active()

    def changed(previous):
        def handler(e):
            if previous:
                previous(e)
            sync()
            page.update()
        return handler

    for control in controls.values():
        if not isinstance(control,(ft.Switch,ft.Checkbox,ft.Dropdown)):
            continue
        # Flet dropdowns emit selection events; switches emit change events.
        event = 'on_select' if hasattr(control, 'on_select') else 'on_change'
        if hasattr(control, event):
            setattr(control, event, changed(getattr(control, event)))
    sync()
    return sync
