type FreeCell<H> = { handle: H | null; free: (h: H) => void };

class HandleRegistry<H> {
  private readonly registry: FinalizationRegistry<FreeCell<H>> | null;

  constructor() {
    this.registry =
      typeof FinalizationRegistry !== "undefined"
        ? new FinalizationRegistry<FreeCell<H>>((cell) => {
            const h = cell.handle;
            if (h == null) return;
            cell.handle = null;
            try {
              cell.free(h);
            } catch (e: unknown) {
              console.warn("Error in HandleRegistry:", e);
            }
          })
        : null;
  }

  register(owner: object, handle: H, free: (h: H) => void): FreeCell<H> {
    const cell: FreeCell<H> = { handle, free };
    this.registry?.register(owner, cell, owner);
    return cell;
  }

  release(owner: object, cell: FreeCell<H>): void {
    const h = cell.handle;
    if (h == null) return;
    cell.free(h);
    cell.handle = null;
    this.registry?.unregister(owner);
  }
}

export const fileHandleRegistry = new HandleRegistry<unknown>();
